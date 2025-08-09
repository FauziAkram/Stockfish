/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2025 The Stockfish developers (see AUTHORS file)

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "search.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <list>
#include <ratio>
#include <string>
#include <utility>

#include "bitboard.h"
#include "evaluate.h"
#include "history.h"
#include "misc.h"
#include "movegen.h"
#include "movepick.h"
#include "nnue/network.h"
#include "nnue/nnue_accumulator.h"
#include "position.h"
#include "syzygy/tbprobe.h"
#include "thread.h"
#include "timeman.h"
#include "tt.h"
#include "uci.h"
#include "ucioption.h"

namespace Stockfish {

namespace TB = Tablebases;

void syzygy_extend_pv(const OptionsMap&            options,
                      const Search::LimitsType&    limits,
                      Stockfish::Position&         pos,
                      Stockfish::Search::RootMove& rootMove,
                      Value&                       v);

using namespace Search;

namespace {

constexpr int SEARCHEDLIST_CAPACITY = 32;
using SearchedList                  = ValueList<Move, SEARCHEDLIST_CAPACITY>;

// FIX: Move forward declarations after SearchedList is defined
template<bool PvNode> void update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus);
template<bool PvNode> void update_quiet_histories(const Position& pos, Stack* ss, Search::Worker& workerThread, Move move, int bonus);
template<bool PvNode> void update_all_stats(const Position& pos, Stack* ss, Search::Worker& workerThread, Move bestMove, Square prevSq, SearchedList& quietsSearched, SearchedList& capturesSearched, Depth depth, Move ttMove);


// (*Scalers):
// The values with Scaler asterisks have proven non-linear scaling.
// They are optimized to time controls of 180 + 1.8 and longer,
// so changing them or adding conditions that are similar requires
// tests at these types of time controls.

template<bool PvNode>
int correction_value(const Worker& w, const Position& pos, const Stack* const ss) {
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))
    const Color us    = pos.side_to_move();
    const auto  m     = (ss - 1)->currentMove;
    const auto  pcv   = w.pawnCorrectionHistory[pawn_correction_history_index(pos)][us];
    const auto  micv  = w.minorPieceCorrectionHistory[minor_piece_index(pos)][us];
    const auto  wnpcv = w.nonPawnCorrectionHistory[non_pawn_index<WHITE>(pos)][WHITE][us];
    const auto  bnpcv = w.nonPawnCorrectionHistory[non_pawn_index<BLACK>(pos)][BLACK][us];
    const auto  cntcv =
      m.is_ok() ? (*(ss - 2)->continuationCorrectionHistory)[pos.piece_on(m.to_sq())][m.to_sq()]
                 : S(0, 0);

    return S(8867, 8867) * pcv + S(8136, 8136) * micv + S(10757, 10757) * (wnpcv + bnpcv) + S(7232, 7232) * cntcv;
#undef S
}

// Add correctionHistory value to raw staticEval and guarantee evaluation
// does not hit the tablebase range.
Value to_corrected_static_eval(const Value v, const int cv) {
    return std::clamp(v + cv / 131072, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);
}

template<bool PvNode>
void update_correction_history(const Position& pos,
                               Stack* const    ss,
                               Search::Worker& workerThread,
                               const int       bonus) {
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))
    const Move  m  = (ss - 1)->currentMove;
    const Color us = pos.side_to_move();

    static constexpr int nonPawnWeight = 165;

    workerThread.pawnCorrectionHistory[pawn_correction_history_index(pos)][us] << bonus;
    workerThread.minorPieceCorrectionHistory[minor_piece_index(pos)][us] << bonus * S(153, 153) / S(128, 128);
    workerThread.nonPawnCorrectionHistory[non_pawn_index<WHITE>(pos)][WHITE][us]
      << bonus * S(nonPawnWeight, nonPawnWeight) / S(128, 128);
    workerThread.nonPawnCorrectionHistory[non_pawn_index<BLACK>(pos)][BLACK][us]
      << bonus * S(nonPawnWeight, nonPawnWeight) / S(128, 128);

    if (m.is_ok())
        (*(ss - 2)->continuationCorrectionHistory)[pos.piece_on(m.to_sq())][m.to_sq()]
          << bonus * S(153, 153) / S(128, 128);
#undef S
}

// Add a small random component to draw evaluations to avoid 3-fold blindness
Value value_draw(size_t nodes) { return VALUE_DRAW - 1 + Value(nodes & 0x2); }
Value value_to_tt(Value v, int ply);
Value value_from_tt(Value v, int ply, int r50c);
void  update_pv(Move* pv, Move move, const Move* childPv);


}  // namespace

Search::Worker::Worker(SharedState&                    sharedState,
                       std::unique_ptr<ISearchManager> sm,
                       size_t                          threadId,
                       NumaReplicatedAccessToken       token) :
    // Unpack the SharedState struct into member variables
    threadIdx(threadId),
    numaAccessToken(token),
    manager(std::move(sm)),
    options(sharedState.options),
    threads(sharedState.threads),
    tt(sharedState.tt),
    networks(sharedState.networks),
    refreshTable(networks[token]) {
    clear();
}

void Search::Worker::ensure_network_replicated() {
    // Access once to force lazy initialization.
    // We do this because we want to avoid initialization during search.
    (void) (networks[numaAccessToken]);
}

void Search::Worker::start_searching() {

    accumulatorStack.reset();

    // Non-main threads go directly to iterative_deepening()
    if (!is_mainthread())
    {
        iterative_deepening();
        return;
    }

    main_manager()->tm.init(limits, rootPos.side_to_move(), rootPos.game_ply(), options,
                            main_manager()->originalTimeAdjust);
    tt.new_search();

    if (rootMoves.empty())
    {
        rootMoves.emplace_back(Move::none());
        main_manager()->updates.onUpdateNoMoves(
          {0, {rootPos.checkers() ? -VALUE_MATE : VALUE_DRAW, rootPos}});
    }
    else
    {
        threads.start_searching();  // start non-main threads
        iterative_deepening();      // main thread start searching
    }

    // When we reach the maximum depth, we can arrive here without a raise of
    // threads.stop. However, if we are pondering or in an infinite search,
    // the UCI protocol states that we shouldn't print the best move before the
    // GUI sends a "stop" or "ponderhit" command. We therefore simply wait here
    // until the GUI sends one of those commands.
    while (!threads.stop && (main_manager()->ponder || limits.infinite))
    {}  // Busy wait for a stop or a ponder reset

    // Stop the threads if not already stopped (also raise the stop if
    // "ponderhit" just reset threads.ponder)
    threads.stop = true;

    // Wait until all threads have finished
    threads.wait_for_search_finished();

    // When playing in 'nodes as time' mode, subtract the searched nodes from
    // the available ones before exiting.
    if (limits.npmsec)
        main_manager()->tm.advance_nodes_time(threads.nodes_searched()
                                              - limits.inc[rootPos.side_to_move()]);

    Worker* bestThread = this;
    Skill   skill =
      Skill(options["Skill Level"], options["UCI_LimitStrength"] ? int(options["UCI_Elo"]) : 0);

    if (int(options["MultiPV"]) == 1 && !limits.depth && !limits.mate && !skill.enabled()
        && rootMoves[0].pv[0] != Move::none())
        bestThread = threads.get_best_thread()->worker.get();

    main_manager()->bestPreviousScore        = bestThread->rootMoves[0].score;
    main_manager()->bestPreviousAverageScore = bestThread->rootMoves[0].averageScore;

    // Send again PV info if we have a new best thread
    if (bestThread != this)
        main_manager()->pv(*bestThread, threads, tt, bestThread->completedDepth);

    std::string ponder;

    if (bestThread->rootMoves[0].pv.size() > 1
        || bestThread->rootMoves[0].extract_ponder_from_tt(tt, rootPos))
        ponder = UCIEngine::move(bestThread->rootMoves[0].pv[1], rootPos.is_chess960());

    auto bestmove = UCIEngine::move(bestThread->rootMoves[0].pv[0], rootPos.is_chess960());
    main_manager()->updates.onBestmove(bestmove, ponder);
}

// Main iterative deepening loop. It calls search()
// repeatedly with increasing depth until the allocated thinking time has been
// consumed, the user stops the search, or the maximum search depth is reached.
void Search::Worker::iterative_deepening() {
#define P(p) (p) // Parameter macro for iterative deepening (not PV-dependent)

    SearchManager* mainThread = (is_mainthread() ? main_manager() : nullptr);

    Move pv[MAX_PLY + 1];

    Depth lastBestMoveDepth = P(0);
    Value lastBestScore     = -VALUE_INFINITE;
    auto  lastBestPV        = std::vector{Move::none()};

    Value  alpha, beta;
    Value  bestValue     = -VALUE_INFINITE;
    Color  us            = rootPos.side_to_move();
    double timeReduction = P(1.0), totBestMoveChanges = P(0.0);
    int    delta, iterIdx = P(0);

    // Allocate stack with extra size to allow access from (ss - 7) to (ss + 2):
    // (ss - 7) is needed for update_continuation_histories(ss - 1) which accesses (ss - 6),
    // (ss + 2) is needed for initialization of cutOffCnt.
    Stack  stack[MAX_PLY + 10] = {};
    Stack* ss                  = stack + 7;

    for (int i = 7; i > 0; --i)
    {
        (ss - i)->continuationHistory =
          &continuationHistory[0][0][NO_PIECE][0];  // Use as a sentinel
        (ss - i)->continuationCorrectionHistory = &continuationCorrectionHistory[NO_PIECE][0];
        (ss - i)->staticEval                    = VALUE_NONE;
    }

    for (int i = 0; i <= MAX_PLY + 2; ++i)
        (ss + i)->ply = i;

    ss->pv = pv;

    if (mainThread)
    {
        if (mainThread->bestPreviousScore == VALUE_INFINITE)
            mainThread->iterValue.fill(VALUE_ZERO);
        else
            mainThread->iterValue.fill(mainThread->bestPreviousScore);
    }

    size_t multiPV = size_t(options["MultiPV"]);
    Skill skill(options["Skill Level"], options["UCI_LimitStrength"] ? int(options["UCI_Elo"]) : 0);

    // When playing with strength handicap enable MultiPV search that we will
    // use behind-the-scenes to retrieve a set of possible moves.
    if (skill.enabled())
        multiPV = std::max(multiPV, size_t(P(4)));

    multiPV = std::min(multiPV, rootMoves.size());

    int searchAgainCounter = P(0);

    lowPlyHistory.fill(P(89));

    // Iterative deepening loop until requested to stop or the target depth is reached
    while (++rootDepth < MAX_PLY && !threads.stop
           && !(limits.depth && mainThread && rootDepth > limits.depth))
    {
        // Age out PV variability metric
        if (mainThread)
            totBestMoveChanges /= P(2.0);

        // Save the last iteration's scores before the first PV line is searched and
        // all the move scores except the (new) PV are set to -VALUE_INFINITE.
        for (RootMove& rm : rootMoves)
            rm.previousScore = rm.score;

        size_t pvFirst = P(0);
        pvLast         = P(0);

        if (!threads.increaseDepth)
            searchAgainCounter++;

        // MultiPV loop. We perform a full root search for each PV line
        for (pvIdx = 0; pvIdx < multiPV; ++pvIdx)
        {
            if (pvIdx == pvLast)
            {
                pvFirst = pvLast;
                for (pvLast++; pvLast < rootMoves.size(); pvLast++)
                    if (rootMoves[pvLast].tbRank != rootMoves[pvFirst].tbRank)
                        break;
            }

            // Reset UCI info selDepth for each depth and each PV line
            selDepth = P(0);

            // Reset aspiration window starting size
            delta     = P(5) + std::abs(rootMoves[pvIdx].meanSquaredScore) / P(11131);
            Value avg = rootMoves[pvIdx].averageScore;
            alpha     = std::max(avg - delta, -VALUE_INFINITE);
            beta      = std::min(avg + delta, VALUE_INFINITE);

            // Adjust optimism based on root move's averageScore
            optimism[us]  = P(136) * avg / (std::abs(avg) + P(93));
            optimism[~us] = -optimism[us];

            // Start with a small aspiration window and, in the case of a fail
            // high/low, re-search with a bigger window until we don't fail
            // high/low anymore.
            int failedHighCnt = P(0);
            while (true)
            {
                // Adjust the effective depth searched, but ensure at least one
                // effective increment for every four searchAgain steps (see issue #2717).
                Depth adjustedDepth =
                  std::max(Depth(P(1)), rootDepth - failedHighCnt - P(3) * (searchAgainCounter + P(1)) / P(4));
                rootDelta = beta - alpha;
                bestValue = search<Root>(rootPos, ss, alpha, beta, adjustedDepth, false);

                // Bring the best move to the front. It is critical that sorting
                // is done with a stable algorithm because all the values but the
                // first and eventually the new best one is set to -VALUE_INFINITE
                // and we want to keep the same order for all the moves except the
                // new PV that goes to the front. Note that in the case of MultiPV
                // search the already searched PV lines are preserved.
                std::stable_sort(rootMoves.begin() + pvIdx, rootMoves.begin() + pvLast);

                // If search has been stopped, we break immediately. Sorting is
                // safe because RootMoves is still valid, although it refers to
                // the previous iteration.
                if (threads.stop)
                    break;

                // When failing high/low give some update before a re-search. To avoid
                // excessive output that could hang GUIs like Fritz 19, only start
                // at nodes > 10M (rather than depth N, which can be reached quickly)
                if (mainThread && multiPV == P(1) && (bestValue <= alpha || bestValue >= beta)
                    && nodes > P(10000000))
                    main_manager()->pv(*this, threads, tt, rootDepth);

                // In case of failing low/high increase aspiration window and re-search,
                // otherwise exit the loop.
                if (bestValue <= alpha)
                {
                    beta  = (P(3) * alpha + beta) / P(4);
                    alpha = std::max(bestValue - delta, -VALUE_INFINITE);

                    failedHighCnt = P(0);
                    if (mainThread)
                        mainThread->stopOnPonderhit = false;
                }
                else if (bestValue >= beta)
                {
                    beta = std::min(bestValue + delta, VALUE_INFINITE);
                    ++failedHighCnt;
                }
                else
                    break;

                delta += delta / P(3);

                assert(alpha >= -VALUE_INFINITE && beta <= VALUE_INFINITE);
            }

            // Sort the PV lines searched so far and update the GUI
            std::stable_sort(rootMoves.begin() + pvFirst, rootMoves.begin() + pvIdx + 1);

            if (mainThread
                && (threads.stop || pvIdx + 1 == multiPV || nodes > P(10000000))
                // A thread that aborted search can have mated-in/TB-loss PV and
                // score that cannot be trusted, i.e. it can be delayed or refuted
                // if we would have had time to fully search other root-moves. Thus
                // we suppress this output and below pick a proven score/PV for this
                // thread (from the previous iteration).
                && !(threads.abortedSearch && is_loss(rootMoves[0].uciScore)))
                main_manager()->pv(*this, threads, tt, rootDepth);

            if (threads.stop)
                break;
        }

        if (!threads.stop)
            completedDepth = rootDepth;

        // We make sure not to pick an unproven mated-in score,
        // in case this thread prematurely stopped search (aborted-search).
        if (threads.abortedSearch && rootMoves[0].score != -VALUE_INFINITE
            && is_loss(rootMoves[0].score))
        {
            // Bring the last best move to the front for best thread selection.
            Utility::move_to_front(rootMoves, [&lastBestPV = std::as_const(lastBestPV)](
                                                const auto& rm) { return rm == lastBestPV[0]; });
            rootMoves[0].pv    = lastBestPV;
            rootMoves[0].score = rootMoves[0].uciScore = lastBestScore;
        }
        else if (rootMoves[0].pv[0] != lastBestPV[0])
        {
            lastBestPV        = rootMoves[0].pv;
            lastBestScore     = rootMoves[0].score;
            lastBestMoveDepth = rootDepth;
        }

        if (!mainThread)
            continue;

        // Have we found a "mate in x"?
        if (limits.mate && rootMoves[0].score == rootMoves[0].uciScore
            && ((rootMoves[0].score >= VALUE_MATE_IN_MAX_PLY
                 && VALUE_MATE - rootMoves[0].score <= P(2) * limits.mate)
                || (rootMoves[0].score != -VALUE_INFINITE
                    && rootMoves[0].score <= VALUE_MATED_IN_MAX_PLY
                    && VALUE_MATE + rootMoves[0].score <= P(2) * limits.mate)))
            threads.stop = true;

        // If the skill level is enabled and time is up, pick a sub-optimal best move
        if (skill.enabled() && skill.time_to_pick(rootDepth))
            skill.pick_best(rootMoves, multiPV);

        // Use part of the gained time from a previous stable move for the current move
        for (auto&& th : threads)
        {
            totBestMoveChanges += th->worker->bestMoveChanges;
            th->worker->bestMoveChanges = 0;
        }

        // Do we have time for the next iteration? Can we stop searching now?
        if (limits.use_time_management() && !threads.stop && !mainThread->stopOnPonderhit)
        {
            uint64_t nodesEffort =
              rootMoves[0].effort * P(100000) / std::max(size_t(P(1)), size_t(nodes));

            double fallingEval =
              (P(11.396) + P(2.035) * (mainThread->bestPreviousAverageScore - bestValue)
               + P(0.968) * (mainThread->iterValue[iterIdx] - bestValue))
              / P(100.0);
            fallingEval = std::clamp(fallingEval, P(0.5786), P(1.6752));

            // If the bestMove is stable over several iterations, reduce time accordingly
            double k      = P(0.527);
            double center = lastBestMoveDepth + P(11.0);
            timeReduction = P(0.8) + P(0.84) / (P(1.077) + std::exp(-k * (completedDepth - center)));
            double reduction =
              (P(1.4540) + mainThread->previousTimeReduction) / (P(2.1593) * timeReduction);
            double bestMoveInstability = P(0.9929) + P(1.8519) * totBestMoveChanges / threads.size();

            double totalTime =
              mainThread->tm.optimum() * fallingEval * reduction * bestMoveInstability;

            // Cap used time in case of a single legal move for a better viewer experience
            if (rootMoves.size() == 1)
                totalTime = std::min(P(500.0), totalTime);

            auto elapsedTime = elapsed();

            if (completedDepth >= P(10) && nodesEffort >= P(97056) && elapsedTime > totalTime * P(0.6540)
                && !mainThread->ponder)
                threads.stop = true;

            // Stop the search if we have exceeded the totalTime or maximum
            if (elapsedTime > std::min(totalTime, double(mainThread->tm.maximum())))
            {
                // If we are allowed to ponder do not stop the search now but
                // keep pondering until the GUI sends "ponderhit" or "stop".
                if (mainThread->ponder)
                    mainThread->stopOnPonderhit = true;
                else
                    threads.stop = true;
            }
            else
                threads.increaseDepth = mainThread->ponder || elapsedTime <= totalTime * P(0.5138);
        }

        mainThread->iterValue[iterIdx] = bestValue;
        iterIdx                        = (iterIdx + 1) & 3;
    }

    if (!mainThread)
        return;

    mainThread->previousTimeReduction = timeReduction;

    // If the skill level is enabled, swap the best PV line with the sub-optimal one
    if (skill.enabled())
        std::swap(rootMoves[0],
                  *std::find(rootMoves.begin(), rootMoves.end(),
                             skill.best ? skill.best : skill.pick_best(rootMoves, multiPV)));
#undef P
}


void Search::Worker::do_move(Position& pos, const Move move, StateInfo& st, Stack* const ss) {
    do_move(pos, move, st, pos.gives_check(move), ss);
}

void Search::Worker::do_move(
  Position& pos, const Move move, StateInfo& st, const bool givesCheck, Stack* const ss) {
    bool       capture = pos.capture_stage(move);
    DirtyPiece dp      = pos.do_move(move, st, givesCheck, &tt);
    nodes.fetch_add(1, std::memory_order_relaxed);
    accumulatorStack.push(dp);
    if (ss != nullptr)
    {
        ss->currentMove         = move;
        ss->continuationHistory = &continuationHistory[ss->inCheck][capture][dp.pc][move.to_sq()];
        ss->continuationCorrectionHistory = &continuationCorrectionHistory[dp.pc][move.to_sq()];
    }
}

void Search::Worker::do_null_move(Position& pos, StateInfo& st) { pos.do_null_move(st, tt); }

void Search::Worker::undo_move(Position& pos, const Move move) {
    pos.undo_move(move);
    accumulatorStack.pop();
}

void Search::Worker::undo_null_move(Position& pos) { pos.undo_null_move(); }


// Reset histories, usually before a new game
void Search::Worker::clear() {
    mainHistory.fill(64);
    captureHistory.fill(-753);
    pawnHistory.fill(-1275);
    pawnCorrectionHistory.fill(5);
    minorPieceCorrectionHistory.fill(0);
    nonPawnCorrectionHistory.fill(0);

    ttMoveHistory = 0;

    for (auto& to : continuationCorrectionHistory)
        for (auto& h : to)
            h.fill(8);

    for (bool inCheck : {false, true})
        for (StatsType c : {NoCaptures, Captures})
            for (auto& to : continuationHistory[inCheck][c])
                for (auto& h : to)
                    h.fill(-494);

    for (size_t i = 1; i < reductions.size(); ++i)
        reductions[i] = int(2782 / 128.0 * std::log(i));

    refreshTable.clear(networks[numaAccessToken]);
}


// Main search function for both PV and non-PV nodes
template<NodeType nodeType>
Value Search::Worker::search(
  Position& pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode) {

    constexpr bool PvNode   = nodeType != NonPV;
    constexpr bool rootNode = nodeType == Root;
    const bool     allNode  = !(PvNode || cutNode);

#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))

    // Dive into quiescence search when the depth reaches zero
    if (depth <= S(0, 0))
    {
        constexpr auto nt = PvNode ? PV : NonPV;
        return qsearch<nt>(pos, ss, alpha, beta);
    }

    // Limit the depth if extensions made it too large
    depth = std::min(depth, Depth(MAX_PLY - S(1, 1)));

    // Check if we have an upcoming move that draws by repetition
    if (!rootNode && alpha < VALUE_DRAW && pos.upcoming_repetition(ss->ply))
    {
        alpha = value_draw(nodes);
        if (alpha >= beta)
            return alpha;
    }

    assert(-VALUE_INFINITE <= alpha && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - S(1, 1)));
    assert(S(0, 0) < depth && depth < MAX_PLY);
    assert(!(PvNode && cutNode));

    Move      pv[MAX_PLY + 1];
    StateInfo st;

    Key   posKey;
    Move  move, excludedMove, bestMove;
    Depth extension, newDepth;
    Value bestValue, value, eval, maxValue, probCutBeta;
    bool  givesCheck, improving, priorCapture, opponentWorsening;
    bool  capture, ttCapture;
    int   priorReduction;
    Piece movedPiece;

    SearchedList capturesSearched;
    SearchedList quietsSearched;

    // Step 1. Initialize node
    ss->inCheck   = pos.checkers();
    priorCapture  = pos.captured_piece();
    Color us      = pos.side_to_move();
    ss->moveCount = S(0, 0);
    bestValue     = -VALUE_INFINITE;
    maxValue      = VALUE_INFINITE;

    // Check for the available remaining time
    if (is_mainthread())
        main_manager()->check_time(*this);

    // Used to send selDepth info to GUI (selDepth counts from 1, ply from 0)
    if (PvNode && selDepth < ss->ply + S(1, 1))
        selDepth = ss->ply + S(1, 1);

    if (!rootNode)
    {
        // Step 2. Check for aborted search and immediate draw
        if (threads.stop.load(std::memory_order_relaxed) || pos.is_draw(ss->ply)
            || ss->ply >= MAX_PLY)
            return (ss->ply >= MAX_PLY && !ss->inCheck) ? evaluate(pos) : value_draw(nodes);

        // Step 3. Mate distance pruning.
        alpha = std::max(mated_in(ss->ply), alpha);
        beta  = std::min(mate_in(ss->ply + S(1, 1)), beta);
        if (alpha >= beta)
            return alpha;
    }

    assert(S(0, 0) <= ss->ply && ss->ply < MAX_PLY);

    Square prevSq  = ((ss - 1)->currentMove).is_ok() ? ((ss - 1)->currentMove).to_sq() : SQ_NONE;
    bestMove       = Move::none();
    priorReduction = (ss - 1)->reduction;
    (ss - 1)->reduction = S(0, 0);
    ss->statScore       = S(0, 0);
    (ss + 2)->cutoffCnt = S(0, 0);

    // Step 4. Transposition table lookup
    excludedMove                   = ss->excludedMove;
    posKey                         = pos.key();
    auto [ttHit, ttData, ttWriter] = tt.probe(posKey);
    // Need further processing of the saved data
    ss->ttHit    = ttHit;
    ttData.move  = rootNode ? rootMoves[pvIdx].pv[0] : ttHit ? ttData.move : Move::none();
    ttData.value = ttHit ? value_from_tt(ttData.value, ss->ply, pos.rule50_count()) : VALUE_NONE;
    ss->ttPv     = excludedMove ? ss->ttPv : PvNode || (ttHit && ttData.is_pv);
    ttCapture    = ttData.move && pos.capture_stage(ttData.move);

    // At non-PV nodes we check for an early TT cutoff
    if (!PvNode && !excludedMove && ttData.depth > depth - S(1, 1) * (ttData.value <= beta)
        && is_valid(ttData.value)
        && (ttData.bound & (ttData.value >= beta ? BOUND_LOWER : BOUND_UPPER))
        && (cutNode == (ttData.value >= beta) || depth > S(5, 5)))
    {
        if (ttData.move && ttData.value >= beta)
        {
            if (!ttCapture)
                update_quiet_histories<PvNode>(pos, ss, *this, ttData.move,
                                       std::min(S(127, 127) * depth - S(74, 74), S(1063, 1063)));

            if (prevSq != SQ_NONE && (ss - 1)->moveCount <= S(3, 3) && !priorCapture)
                update_continuation_histories<PvNode>(ss - 1, pos.piece_on(prevSq), prevSq, S(-2128, -2128));
        }

        if (pos.rule50_count() < S(91, 91))
        {
            if (depth >= S(8, 8) && ttData.move && pos.pseudo_legal(ttData.move) && pos.legal(ttData.move)
                && !is_decisive(ttData.value))
            {
                pos.do_move(ttData.move, st);
                Key nextPosKey                             = pos.key();
                auto [ttHitNext, ttDataNext, ttWriterNext] = tt.probe(nextPosKey);
                pos.undo_move(ttData.move);

                if (!is_valid(ttDataNext.value))
                    return ttData.value;
                if ((ttData.value >= beta) == (-ttDataNext.value >= beta))
                    return ttData.value;
            }
            else
                return ttData.value;
        }
    }

    // Step 5. Tablebases probe
    if (!rootNode && !excludedMove && tbConfig.cardinality)
    {
        int piecesCount = pos.count<ALL_PIECES>();

        if (piecesCount <= tbConfig.cardinality
            && (piecesCount < tbConfig.cardinality || depth >= tbConfig.probeDepth)
            && pos.rule50_count() == S(0, 0) && !pos.can_castle(ANY_CASTLING))
        {
            TB::ProbeState err;
            TB::WDLScore   wdl = Tablebases::probe_wdl(pos, &err);

            if (is_mainthread())
                main_manager()->callsCnt = S(0, 0);

            if (err != TB::ProbeState::FAIL)
            {
                tbHits.fetch_add(S(1, 1), std::memory_order_relaxed);

                int drawScore = tbConfig.useRule50 ? S(1, 1) : S(0, 0);

                Value tbValue = VALUE_TB - ss->ply;

                value = wdl < -drawScore ? -tbValue
                      : wdl > drawScore  ? tbValue
                                         : VALUE_DRAW + S(2, 2) * wdl * drawScore;
                Bound b = wdl < -drawScore ? BOUND_UPPER
                        : wdl > drawScore  ? BOUND_LOWER
                                           : BOUND_EXACT;

                if (b == BOUND_EXACT || (b == BOUND_LOWER ? value >= beta : value <= alpha))
                {
                    ttWriter.write(posKey, value_to_tt(value, ss->ply), ss->ttPv, b,
                                   std::min(Depth(MAX_PLY - S(1, 1)), depth + S(6, 6)), Move::none(), VALUE_NONE,
                                   tt.generation());

                    return value;
                }
                if (PvNode)
                {
                    if (b == BOUND_LOWER)
                        bestValue = value, alpha = std::max(alpha, bestValue);
                    else
                        maxValue = value;
                }
            }
        }
    }

    // Step 6. Static evaluation of the position
    Value unadjustedStaticEval = VALUE_NONE;
    const auto correctionValue = correction_value<PvNode>(*this, pos, ss);

    if (ss->inCheck)
    {
        ss->staticEval = eval = (ss - 2)->staticEval;
        improving             = false;
        goto moves_loop;
    }
    else if (excludedMove)
        unadjustedStaticEval = eval = ss->staticEval;
    else if (ss->ttHit)
    {
        unadjustedStaticEval = ttData.eval;
        if (!is_valid(unadjustedStaticEval))
            unadjustedStaticEval = evaluate(pos);
        ss->staticEval = eval = to_corrected_static_eval(unadjustedStaticEval, correctionValue);

        if (is_valid(ttData.value)
            && (ttData.bound & (ttData.value > eval ? BOUND_LOWER : BOUND_UPPER)))
            eval = ttData.value;
    }
    else
    {
        unadjustedStaticEval = evaluate(pos);
        ss->staticEval = eval = to_corrected_static_eval(unadjustedStaticEval, correctionValue);
        ttWriter.write(posKey, VALUE_NONE, ss->ttPv, BOUND_NONE, DEPTH_UNSEARCHED, Move::none(),
                       unadjustedStaticEval, tt.generation());
    }

    if (((ss - 1)->currentMove).is_ok() && !(ss - 1)->inCheck && !priorCapture)
    {
        int bonus = std::clamp(-S(10, 10) * int((ss - 1)->staticEval + ss->staticEval), S(-1979, -1979), S(1561, 1561)) + S(630, 630);
        mainHistory[~us][((ss - 1)->currentMove).from_to()] << bonus * S(935, 935) / S(1024, 1024);
        if (!ttHit && type_of(pos.piece_on(prevSq)) != PAWN
            && ((ss - 1)->currentMove).type_of() != PROMOTION)
            pawnHistory[pawn_history_index(pos)][pos.piece_on(prevSq)][prevSq]
              << bonus * S(1428, 1428) / S(1024, 1024);
    }

    improving         = ss->staticEval > (ss - 2)->staticEval;
    opponentWorsening = ss->staticEval > -(ss - 1)->staticEval;

    if (priorReduction >= (depth < S(10, 10) ? S(1, 1) : S(3, 3)) && !opponentWorsening)
        depth++;
    if (priorReduction >= S(2, 2) && depth >= S(2, 2) && ss->staticEval + (ss - 1)->staticEval > S(177, 177))
        depth--;

    // Step 7. Razoring
    if (!PvNode && eval < alpha - S(495, 495) - S(290, 290) * depth * depth)
        return qsearch<NonPV>(pos, ss, alpha, beta);

    // Step 8. Futility pruning: child node
    {
        auto futility_margin = [&](Depth d) {
            Value futilityMult = S(90, 90) - S(20, 20) * (cutNode && !ss->ttHit);
            return futilityMult * d
                 - improving * futilityMult * S(2, 2)
                 - opponentWorsening * futilityMult / S(3, 3)
                 + (ss - 1)->statScore / S(356, 356)
                 + std::abs(correctionValue) / S(171290, 171290);
        };
        if (!ss->ttPv && depth < S(14, 14) && eval - futility_margin(depth) >= beta && eval >= beta
            && (!ttData.move || ttCapture) && !is_loss(beta) && !is_win(eval))
            return beta + (eval - beta) / S(3, 3);
    }

    // Step 9. Null move search with verification search
    if (cutNode && ss->staticEval >= beta - S(19, 19) * depth + S(403, 403) && !excludedMove
        && pos.non_pawn_material(us) && ss->ply >= nmpMinPly && !is_loss(beta))
    {
        assert((ss - 1)->currentMove != Move::null());

        Depth R = S(7, 7) + depth / S(3, 3);
        ss->currentMove                   = Move::null();
        ss->continuationHistory           = &continuationHistory[0][0][NO_PIECE][0];
        ss->continuationCorrectionHistory = &continuationCorrectionHistory[NO_PIECE][0];
        do_null_move(pos, st);
        Value nullValue = -search<NonPV>(pos, ss + 1, -beta, -beta + S(1, 1), depth - R, false);
        undo_null_move(pos);

        if (nullValue >= beta && !is_win(nullValue))
        {
            if (nmpMinPly || depth < S(16, 16))
                return nullValue;

            assert(!nmpMinPly);
            nmpMinPly = ss->ply + S(3, 3) * (depth - R) / S(4, 4);
            Value v = search<NonPV>(pos, ss, beta - S(1, 1), beta, depth - R, false);
            nmpMinPly = S(0, 0);
            if (v >= beta)
                return nullValue;
        }
    }

    improving |= ss->staticEval >= beta;

    // Step 10. Internal iterative reductions
    if (!allNode && depth >= S(6, 6) && !ttData.move && priorReduction <= S(3, 3))
        depth--;

    // Step 11. ProbCut
    probCutBeta = beta + S(215, 215) - S(60, 60) * improving;
    if (depth >= S(3, 3) && !is_decisive(beta)
        && !(is_valid(ttData.value) && ttData.value < probCutBeta))
    {
        assert(probCutBeta < VALUE_INFINITE && probCutBeta > beta);
        MovePicker mp(pos, ttData.move, probCutBeta - ss->staticEval, &captureHistory);
        Depth      dynamicReduction = (ss->staticEval - beta) / S(300, 300);
        Depth      probCutDepth     = std::max(Depth(S(0, 0)), depth - S(5, 5) - dynamicReduction);

        while ((move = mp.next_move()) != Move::none())
        {
            assert(move.is_ok());
            if (move == excludedMove || !pos.legal(move))
                continue;

            assert(pos.capture_stage(move));
            movedPiece = pos.moved_piece(move);
            do_move(pos, move, st, ss);
            value = -qsearch<NonPV>(pos, ss + 1, -probCutBeta, -probCutBeta + S(1, 1));

            if (value >= probCutBeta && probCutDepth > S(0, 0))
                value = -search<NonPV>(pos, ss + 1, -probCutBeta, -probCutBeta + S(1, 1), probCutDepth, !cutNode);

            undo_move(pos, move);
            if (value >= probCutBeta)
            {
                ttWriter.write(posKey, value_to_tt(value, ss->ply), ss->ttPv, BOUND_LOWER,
                               probCutDepth + S(1, 1), move, unadjustedStaticEval, tt.generation());
                if (!is_decisive(value))
                    return value - (probCutBeta - beta);
            }
        }
    }

moves_loop:
    // Step 12. A small Probcut idea
    probCutBeta = beta + S(417, 417);
    if ((ttData.bound & BOUND_LOWER) && ttData.depth >= depth - S(4, 4) && ttData.value >= probCutBeta
        && !is_decisive(beta) && is_valid(ttData.value) && !is_decisive(ttData.value))
        return probCutBeta;

    const PieceToHistory* contHist[] = {
      (ss - 1)->continuationHistory, (ss - 2)->continuationHistory, (ss - 3)->continuationHistory,
      (ss - 4)->continuationHistory, (ss - 5)->continuationHistory, (ss - 6)->continuationHistory};


    MovePicker mp(pos, ttData.move, depth, &mainHistory, &lowPlyHistory, &captureHistory, contHist,
                  &pawnHistory, ss->ply);

    value = bestValue;
    int moveCount = S(0, 0);

    // Step 13. Loop through all pseudo-legal moves
    while ((move = mp.next_move()) != Move::none())
    {
        assert(move.is_ok());
        if (move == excludedMove) continue;
        if (!pos.legal(move)) continue;
        if (rootNode && !std::count(rootMoves.begin() + pvIdx, rootMoves.begin() + pvLast, move)) continue;

        ss->moveCount = ++moveCount;

        if (rootNode && is_mainthread() && nodes > S(10000000, 10000000))
        {
            main_manager()->updates.onIter(
              {depth, UCIEngine::move(move, pos.is_chess960()), moveCount + pvIdx});
        }
        if (PvNode)
            (ss + 1)->pv = nullptr;

        extension  = S(0, 0);
        capture    = pos.capture_stage(move);
        movedPiece = pos.moved_piece(move);
        givesCheck = pos.gives_check(move);
        (ss + 1)->quietMoveStreak = (!capture && !givesCheck) ? (ss->quietMoveStreak + S(1, 1)) : S(0, 0);

        newDepth = depth - S(1, 1);
        int delta = beta - alpha;
        Depth r = this->reduction<PvNode>(improving, depth, moveCount, delta);

        if (ss->ttPv)
            r += S(931, 931);

        // Step 14. Pruning at shallow depth.
        if (!rootNode && pos.non_pawn_material(us) && !is_loss(bestValue))
        {
            if (moveCount >= (S(3, 3) + depth * depth) / (S(2, 2) - improving))
                mp.skip_quiet_moves();

            int lmrDepth = newDepth - r / S(1024, 1024);

            if (capture || givesCheck)
            {
                Piece capturedPiece = pos.piece_on(move.to_sq());
                int   captHist = captureHistory[movedPiece][move.to_sq()][type_of(capturedPiece)];

                if (!givesCheck && lmrDepth < S(7, 7) && !ss->inCheck)
                {
                    Value futilityValue = ss->staticEval + S(232, 232) + S(224, 224) * lmrDepth
                                        + PieceValue[capturedPiece] + S(131, 131) * captHist / S(1024, 1024);
                    if (futilityValue <= alpha)
                        continue;
                }

                int margin = std::clamp(S(158, 158) * depth + captHist / S(31, 31), S(0, 0), S(283, 283) * depth);
                if (!pos.see_ge(move, -margin))
                {
                    bool mayStalemateTrap =
                      depth > S(2, 2) && alpha < S(0, 0) && pos.non_pawn_material(us) == PieceValue[movedPiece]
                      && PieceValue[movedPiece] >= RookValue
                      && !(attacks_bb<KING>(pos.square<KING>(us)) & move.from_sq())
                      && !mp.can_move_king_or_pawn();
                    if (!mayStalemateTrap)
                        continue;
                }
            }
            else
            {
                int history = (*contHist[0])[movedPiece][move.to_sq()]
                            + (*contHist[1])[movedPiece][move.to_sq()]
                            + pawnHistory[pawn_history_index(pos)][movedPiece][move.to_sq()];

                if (history < -S(4361, 4361) * depth)
                    continue;

                history += S(71, 71) * mainHistory[us][move.from_to()] / S(32, 32);
                lmrDepth += history / S(3233, 3233);
                Value baseFutility = (bestMove ? S(46, 46) : S(230, 230));
                Value futilityValue =
                  ss->staticEval + baseFutility + S(131, 131) * lmrDepth + S(91, 91) * (ss->staticEval > alpha);

                if (!ss->inCheck && lmrDepth < S(11, 11) && futilityValue <= alpha)
                {
                    if (bestValue <= futilityValue && !is_decisive(bestValue)
                        && !is_win(futilityValue))
                        bestValue = futilityValue;
                    continue;
                }

                lmrDepth = std::max(lmrDepth, S(0, 0));
                if (!pos.see_ge(move, -S(26, 26) * lmrDepth * lmrDepth))
                    continue;
            }
        }

        // Step 15. Singular extension search.
        if (!rootNode && move == ttData.move && !excludedMove
            && depth >= S(6, 6) - (completedDepth > S(26, 26)) + ss->ttPv && is_valid(ttData.value)
            && !is_decisive(ttData.value) && (ttData.bound & BOUND_LOWER)
            && ttData.depth >= depth - S(3, 3))
        {
            Value singularBeta  = ttData.value - (S(56, 56) + S(79, 79) * (ss->ttPv && !PvNode)) * depth / S(58, 58);
            Depth singularDepth = newDepth / S(2, 2);

            ss->excludedMove = move;
            value = search<NonPV>(pos, ss, singularBeta - S(1, 1), singularBeta, singularDepth, cutNode);
            ss->excludedMove = Move::none();

            if (value < singularBeta)
            {
                int corrValAdj   = std::abs(correctionValue) / S(249096, 249096);
                int doubleMargin = S(4, 209) - S(223, 223) * !ttCapture - corrValAdj
                                 - S(959, 959) * ttMoveHistory / S(131072, 131072) - (ss->ply > rootDepth) * S(45, 45);
                int tripleMargin = S(80, 356) - S(249, 249) * !ttCapture + S(86, 86) * ss->ttPv - corrValAdj
                                 - (ss->ply * S(2, 2) > rootDepth * S(3, 3)) * S(53, 53);
                extension =
                  S(1, 1) + (value < singularBeta - doubleMargin) + (value < singularBeta - tripleMargin);

                depth++;
            }
            else if (value >= beta && !is_decisive(value))
                return value;
            else if (ttData.value >= beta)
                extension = S(-3, -3);
            else if (cutNode)
                extension = S(-2, -2);
        }

        // Step 16. Make the move
        do_move(pos, move, st, givesCheck, ss);

        newDepth += extension;
        uint64_t nodeCount = rootNode ? uint64_t(nodes) : S(0, 0);

        if (ss->ttPv)
            r -= S(2510, 3473) + (ttData.value > alpha) * S(916, 916)
               + (ttData.depth >= depth) * (S(943, 943) + cutNode * S(1180, 1180));
        
        r += S(650, 650);
        r -= moveCount * S(69, 69);
        r -= std::abs(correctionValue) / S(27160, 27160);
        if (cutNode) r += S(3000, 3000) + S(1024, 1024) * !ttData.move;
        if (ttCapture) r += S(1350, 1350);
        if ((ss + 1)->cutoffCnt > S(2, 2)) r += S(935, 935) + allNode * S(763, 763);
        r += (ss + 1)->quietMoveStreak * S(51, 51);
        if (move == ttData.move) r -= S(2043, 2043);

        if (capture)
            ss->statScore = S(782, 782) * int(PieceValue[pos.captured_piece()]) / S(128, 128)
                          + captureHistory[movedPiece][move.to_sq()][type_of(pos.captured_piece())];
        else
            ss->statScore = S(2, 2) * mainHistory[us][move.from_to()]
                          + (*contHist[0])[movedPiece][move.to_sq()]
                          + (*contHist[1])[movedPiece][move.to_sq()];
        r -= ss->statScore * S(789, 789) / S(8192, 8192);

        // Step 17. Late moves reduction / extension (LMR)
        if (depth >= S(2, 2) && moveCount > S(1, 1))
        {
            Depth d = std::max(Depth(S(1, 1)), std::min(newDepth - r / S(1024, 1024), newDepth + S(1, 2))) + S(0, 1);
            ss->reduction = newDepth - d;
            value         = -search<NonPV>(pos, ss + 1, -(alpha + S(1, 1)), -alpha, d, true);
            ss->reduction = S(0, 0);

            if (value > alpha)
            {
                const bool doDeeperSearch = d < newDepth && value > (bestValue + S(43, 43) + S(2, 2) * newDepth);
                const bool doShallowerSearch = value < bestValue + S(9, 9);
                newDepth += doDeeperSearch - doShallowerSearch;
                if (newDepth > d)
                    value = -search<NonPV>(pos, ss + 1, -(alpha + S(1, 1)), -alpha, newDepth, !cutNode);
                update_continuation_histories<PvNode>(ss, movedPiece, move.to_sq(), S(1412, 1412));
            }
        }
        // Step 18. Full-depth search when LMR is skipped
        else if (!PvNode || moveCount > S(1, 1))
        {
            if (!ttData.move) r += S(1139, 1139);
            const int threshold1 = depth <= S(4, 4) ? S(2000, 2000) : S(3200, 3200);
            const int threshold2 = depth <= S(4, 4) ? S(3500, 3500) : S(4600, 4600);
            value = -search<NonPV>(pos, ss + 1, -(alpha + S(1, 1)), -alpha,
                                   newDepth - (r > threshold1) - (r > threshold2 && newDepth > S(2, 2)), !cutNode);
        }

        // For PV nodes only, do a full PV search
        if (PvNode && (moveCount == S(1, 1) || value > alpha))
        {
            (ss + 1)->pv    = pv;
            (ss + 1)->pv[0] = Move::none();
            if (move == ttData.move && rootDepth > S(8, 8)) newDepth = std::max(newDepth, Depth(S(1, 1)));
            value = -search<PV>(pos, ss + 1, -beta, -alpha, newDepth, false);
        }

        // Step 19. Undo move
        undo_move(pos, move);
        assert(value > -VALUE_INFINITE && value < VALUE_INFINITE);

        // Step 20. Check for a new best move
        if (threads.stop.load(std::memory_order_relaxed)) return VALUE_ZERO;

        if (rootNode)
        {
            RootMove& rm = *std::find(rootMoves.begin(), rootMoves.end(), move);
            rm.effort += nodes - nodeCount;
            rm.averageScore = rm.averageScore != -VALUE_INFINITE ? (value + rm.averageScore) / S(2, 2) : value;
            rm.meanSquaredScore = rm.meanSquaredScore != -VALUE_INFINITE * VALUE_INFINITE
                                  ? (value * std::abs(value) + rm.meanSquaredScore) / S(2, 2)
                                  : value * std::abs(value);

            if (moveCount == S(1, 1) || value > alpha)
            {
                rm.score = rm.uciScore = value;
                rm.selDepth            = selDepth;
                rm.scoreLowerbound = rm.scoreUpperbound = false;
                if (value >= beta) { rm.scoreLowerbound = true; rm.uciScore = beta; }
                else if (value <= alpha) { rm.scoreUpperbound = true; rm.uciScore = alpha; }
                rm.pv.resize(S(1, 1));
                assert((ss + 1)->pv);
                for (Move* m = (ss + 1)->pv; *m != Move::none(); ++m) rm.pv.push_back(*m);
                if (moveCount > S(1, 1) && !pvIdx) ++bestMoveChanges;
            } else rm.score = -VALUE_INFINITE;
        }

        int inc = (value == bestValue && ss->ply + S(2, 2) >= rootDepth && (int(nodes) & S(14, 14)) == S(0, 0)
                   && !is_win(std::abs(value) + S(1, 1)));

        if (value + inc > bestValue)
        {
            bestValue = value;
            if (value + inc > alpha)
            {
                bestMove = move;
                if (PvNode && !rootNode) update_pv(ss->pv, move, (ss + 1)->pv);
                if (value >= beta)
                {
                    ss->cutoffCnt += (extension < S(2, 2)) || PvNode;
                    assert(value >= beta);
                    break;
                }
                if (depth > S(2, 2) && depth < S(16, 16) && !is_decisive(value)) depth -= S(2, 2);
                assert(depth > S(0, 0));
                alpha = value;
            }
        }
        if (move != bestMove && moveCount <= SEARCHEDLIST_CAPACITY)
        {
            if (capture) capturesSearched.push_back(move);
            else quietsSearched.push_back(move);
        }
    }

    // Step 21. Check for mate and stalemate
    assert(moveCount || !ss->inCheck || excludedMove || !MoveList<LEGAL>(pos).size());

    if (bestValue >= beta && !is_decisive(bestValue) && !is_decisive(alpha))
        bestValue = (bestValue * depth + beta) / (depth + S(1, 1));

    if (!moveCount)
        bestValue = excludedMove ? alpha : ss->inCheck ? mated_in(ss->ply) : VALUE_DRAW;
    else if (bestMove)
    {
        update_all_stats<PvNode>(pos, ss, *this, bestMove, prevSq, quietsSearched, capturesSearched, depth, ttData.move);
        if (!PvNode)
            ttMoveHistory << (bestMove == ttData.move ? S(811, 811) : S(-848, -848));
    }
    else if (!priorCapture && prevSq != SQ_NONE)
    {
        int bonusScale = S(-215, -215);
        bonusScale += std::min(-(ss - 1)->statScore / S(103, 103), S(337, 337));
        bonusScale += std::min(S(64, 64) * depth, S(552, 552));
        bonusScale += S(177, 177) * ((ss - 1)->moveCount > S(8, 8));
        bonusScale += S(141, 141) * (!ss->inCheck && bestValue <= ss->staticEval - S(94, 94));
        bonusScale += S(141, 141) * (!(ss - 1)->inCheck && bestValue <= -(ss - 1)->staticEval - S(76, 76));
        bonusScale = std::max(bonusScale, S(0, 0));
        const int scaledBonus = std::min(S(155, 155) * depth - S(88, 88), S(1416, 1416)) * bonusScale;
        update_continuation_histories<PvNode>(ss - 1, pos.piece_on(prevSq), prevSq, scaledBonus * S(397, 397) / S(32768, 32768));
        mainHistory[~us][((ss - 1)->currentMove).from_to()] << scaledBonus * S(224, 224) / S(32768, 32768);
        if (type_of(pos.piece_on(prevSq)) != PAWN && ((ss - 1)->currentMove).type_of() != PROMOTION)
            pawnHistory[pawn_history_index(pos)][pos.piece_on(prevSq)][prevSq]
              << scaledBonus * S(1127, 1127) / S(32768, 32768);
    }
    else if (priorCapture && prevSq != SQ_NONE)
    {
        Piece capturedPiece = pos.captured_piece();
        assert(capturedPiece != NO_PIECE);
        captureHistory[pos.piece_on(prevSq)][prevSq][type_of(capturedPiece)] << S(1042, 1042);
    }

    if (PvNode)
        bestValue = std::min(bestValue, maxValue);

    if (bestValue <= alpha)
        ss->ttPv = ss->ttPv || (ss - 1)->ttPv;

    if (!excludedMove && !(rootNode && pvIdx))
        ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), ss->ttPv,
                       bestValue >= beta    ? BOUND_LOWER
                       : PvNode && bestMove ? BOUND_EXACT
                                            : BOUND_UPPER,
                       moveCount != 0 ? depth : std::min(Depth(MAX_PLY - S(1, 1)), depth + S(6, 6)), bestMove,
                       unadjustedStaticEval, tt.generation());

    if (!ss->inCheck && !(bestMove && pos.capture(bestMove))
        && ((bestValue < ss->staticEval && bestValue < beta)
            || (bestValue > ss->staticEval && bestMove)))
    {
        auto bonus = std::clamp(int(bestValue - ss->staticEval) * depth / S(8, 8),
                                -CORRECTION_HISTORY_LIMIT / S(4, 4), CORRECTION_HISTORY_LIMIT / S(4, 4));
        update_correction_history<PvNode>(pos, ss, *this, bonus);
    }

    assert(bestValue > -VALUE_INFINITE && bestValue < VALUE_INFINITE);
#undef S
    return bestValue;
}


// Quiescence search function
template<NodeType nodeType>
Value Search::Worker::qsearch(Position& pos, Stack* ss, Value alpha, Value beta) {

    static_assert(nodeType != Root);
    constexpr bool PvNode = nodeType == PV;
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))

    assert(alpha >= -VALUE_INFINITE && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - S(1, 1)));

    if (alpha < VALUE_DRAW && pos.upcoming_repetition(ss->ply))
    {
        alpha = value_draw(nodes);
        if (alpha >= beta)
            return alpha;
    }

    Move pv[MAX_PLY + 1];
    StateInfo st;
    Key posKey;
    Move move, bestMove;
    Value bestValue, value, futilityBase;
    bool pvHit, givesCheck, capture;
    int moveCount;

    if (PvNode)
    {
        (ss + 1)->pv = pv;
        ss->pv[0]    = Move::none();
    }

    bestMove    = Move::none();
    ss->inCheck = pos.checkers();
    moveCount   = S(0, 0);

    if (PvNode && selDepth < ss->ply + S(1, 1))
        selDepth = ss->ply + S(1, 1);

    if (pos.is_draw(ss->ply) || ss->ply >= MAX_PLY)
        return (ss->ply >= MAX_PLY && !ss->inCheck) ? evaluate(pos) : VALUE_DRAW;

    assert(S(0, 0) <= ss->ply && ss->ply < MAX_PLY);

    posKey                         = pos.key();
    auto [ttHit, ttData, ttWriter] = tt.probe(posKey);
    ss->ttHit    = ttHit;
    ttData.move  = ttHit ? ttData.move : Move::none();
    ttData.value = ttHit ? value_from_tt(ttData.value, ss->ply, pos.rule50_count()) : VALUE_NONE;
    pvHit        = ttHit && ttData.is_pv;

    if (!PvNode && ttData.depth >= DEPTH_QS && is_valid(ttData.value)
        && (ttData.bound & (ttData.value >= beta ? BOUND_LOWER : BOUND_UPPER)))
        return ttData.value;

    Value unadjustedStaticEval = VALUE_NONE;
    if (ss->inCheck)
        bestValue = futilityBase = -VALUE_INFINITE;
    else
    {
        const auto correctionValue = correction_value<PvNode>(*this, pos, ss);
        if (ss->ttHit)
        {
            unadjustedStaticEval = ttData.eval;
            if (!is_valid(unadjustedStaticEval))
                unadjustedStaticEval = evaluate(pos);
            ss->staticEval = bestValue = to_corrected_static_eval(unadjustedStaticEval, correctionValue);
            if (is_valid(ttData.value) && !is_decisive(ttData.value)
                && (ttData.bound & (ttData.value > bestValue ? BOUND_LOWER : BOUND_UPPER)))
                bestValue = ttData.value;
        }
        else
        {
            unadjustedStaticEval = evaluate(pos);
            ss->staticEval = bestValue = to_corrected_static_eval(unadjustedStaticEval, correctionValue);
        }

        if (bestValue >= beta)
        {
            if (!is_decisive(bestValue))
                bestValue = (bestValue + beta) / S(2, 2);
            if (!ss->ttHit)
                ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), false, BOUND_LOWER,
                               DEPTH_UNSEARCHED, Move::none(), unadjustedStaticEval,
                               tt.generation());
            return bestValue;
        }
        if (bestValue > alpha)
            alpha = bestValue;
        futilityBase = ss->staticEval + S(352, 352);
    }

    const PieceToHistory* contHist[] = {(ss - 1)->continuationHistory, (ss - 2)->continuationHistory};
    Square prevSq = ((ss - 1)->currentMove).is_ok() ? ((ss - 1)->currentMove).to_sq() : SQ_NONE;

    MovePicker mp(pos, ttData.move, DEPTH_QS, &mainHistory, &lowPlyHistory, &captureHistory,
                  contHist, &pawnHistory, ss->ply);

    while ((move = mp.next_move()) != Move::none())
    {
        assert(move.is_ok());
        if (!pos.legal(move)) continue;

        givesCheck = pos.gives_check(move);
        capture    = pos.capture_stage(move);
        moveCount++;

        if (!is_loss(bestValue))
        {
            if (!givesCheck && move.to_sq() != prevSq && !is_loss(futilityBase) && move.type_of() != PROMOTION)
            {
                if (moveCount > S(2, 2)) continue;
                Value futilityValue = futilityBase + PieceValue[pos.piece_on(move.to_sq())];
                if (futilityValue <= alpha) { bestValue = std::max(bestValue, futilityValue); continue; }
                if (!pos.see_ge(move, alpha - futilityBase)) { bestValue = std::min(alpha, futilityBase); continue; }
            }
            if (!capture && (*contHist[0])[pos.moved_piece(move)][move.to_sq()]
                   + pawnHistory[pawn_history_index(pos)][pos.moved_piece(move)][move.to_sq()]
                 <= S(5868, 5868))
                continue;
            if (!pos.see_ge(move, S(-74, -74))) continue;
        }
        do_move(pos, move, st, givesCheck, ss);
        value = -qsearch<nodeType>(pos, ss + 1, -beta, -alpha);
        undo_move(pos, move);
        assert(value > -VALUE_INFINITE && value < VALUE_INFINITE);

        if (value > bestValue)
        {
            bestValue = value;
            if (value > alpha)
            {
                bestMove = move;
                if (PvNode) update_pv(ss->pv, move, (ss + 1)->pv);
                if (value < beta) alpha = value;
                else break;
            }
        }
    }

    if (ss->inCheck && bestValue == -VALUE_INFINITE)
    {
        assert(!MoveList<LEGAL>(pos).size());
        return mated_in(ss->ply);
    }
    if (!is_decisive(bestValue) && bestValue > beta)
        bestValue = (bestValue + beta) / S(2, 2);

    Color us = pos.side_to_move();
    if (!ss->inCheck && !moveCount && !pos.non_pawn_material(us) && type_of(pos.captured_piece()) >= ROOK)
    {
        if (!((us == WHITE ? shift<NORTH>(pos.pieces(us, PAWN)) : shift<SOUTH>(pos.pieces(us, PAWN))) & ~pos.pieces()))
        {
            pos.state()->checkersBB = Rank1BB;
            if (!MoveList<LEGAL>(pos).size()) bestValue = VALUE_DRAW;
            pos.state()->checkersBB = S(0, 0);
        }
    }

    ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), pvHit,
                   bestValue >= beta ? BOUND_LOWER : BOUND_UPPER, DEPTH_QS, bestMove,
                   unadjustedStaticEval, tt.generation());

    assert(bestValue > -VALUE_INFINITE && bestValue < VALUE_INFINITE);
#undef S
    return bestValue;
}

// FIX: Define the templated member function reduction() correctly
template<bool PvNode>
Depth Search::Worker::reduction(bool i, Depth d, int mn, int delta) const {
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))
    int reductionScale = reductions[d] * reductions[mn];
    return reductionScale - delta * S(731, 731) / rootDelta + !i * reductionScale * S(216, 216) / S(512, 512) + S(1089, 1089);
#undef S
}

TimePoint Search::Worker::elapsed() const {
    return main_manager()->tm.elapsed([this]() { return threads.nodes_searched(); });
}

TimePoint Search::Worker::elapsed_time() const { return main_manager()->tm.elapsed_time(); }

Value Search::Worker::evaluate(const Position& pos) {
    return Eval::evaluate(networks[numaAccessToken], pos, accumulatorStack, refreshTable,
                          optimism[pos.side_to_move()]);
}

namespace {
// Adjusts a mate or TB score from "plies to mate from the root" to
// "plies to mate from the current position". Standard scores are unchanged.
// The function is called before storing a value in the transposition table.
Value value_to_tt(Value v, int ply) { return is_win(v) ? v + ply : is_loss(v) ? v - ply : v; }


// Inverse of value_to_tt()
Value value_from_tt(Value v, int ply, int r50c) {
    if (!is_valid(v))
        return VALUE_NONE;
    if (is_win(v))
    {
        if (v >= VALUE_MATE_IN_MAX_PLY && VALUE_MATE - v > 100 - r50c) return VALUE_TB_WIN_IN_MAX_PLY - 1;
        if (VALUE_TB - v > 100 - r50c) return VALUE_TB_WIN_IN_MAX_PLY - 1;
        return v - ply;
    }
    if (is_loss(v))
    {
        if (v <= VALUE_MATED_IN_MAX_PLY && VALUE_MATE + v > 100 - r50c) return VALUE_TB_LOSS_IN_MAX_PLY + 1;
        if (VALUE_TB + v > 100 - r50c) return VALUE_TB_LOSS_IN_MAX_PLY + 1;
        return v + ply;
    }
    return v;
}


// Adds current move and appends child pv[]
void update_pv(Move* pv, Move move, const Move* childPv) {
    for (*pv++ = move; childPv && *childPv != Move::none();)
        *pv++ = *childPv++;
    *pv = Move::none();
}


// Updates stats at the end of search() when a bestMove is found
template<bool PvNode>
void update_all_stats(const Position& pos,
                      Stack*          ss,
                      Search::Worker& workerThread,
                      Move            bestMove,
                      Square          prevSq,
                      SearchedList&   quietsSearched,
                      SearchedList&   capturesSearched,
                      Depth           depth,
                      Move            ttMove) {
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))
    CapturePieceToHistory& captureHistory = workerThread.captureHistory;
    Piece                  movedPiece     = pos.moved_piece(bestMove);
    PieceType              capturedPiece;

    int quietBonus   = std::min(S(170, 170) * depth - S(87, 87), S(1598, 1598)) + S(332, 332) * (bestMove == ttMove);
    int quietMalus   = std::min(S(743, 743) * depth - S(180, 180), S(2287, 2287)) - S(33, 33) * quietsSearched.size();
    int captureBonus = std::min(S(124, 124) * depth - S(62, 62), S(1245, 1245)) + S(336, 336) * (bestMove == ttMove);
    int captureMalus = std::min(S(708, 708) * depth - S(148, 148), S(2287, 2287)) - S(29, 29) * capturesSearched.size();

    if (!pos.capture_stage(bestMove))
    {
        update_quiet_histories<PvNode>(pos, ss, workerThread, bestMove, quietBonus * S(978, 978) / S(1024, 1024));
        for (Move move : quietsSearched)
            update_quiet_histories<PvNode>(pos, ss, workerThread, move, -quietMalus * S(1115, 1115) / S(1024, 1024));
    }
    else
    {
        capturedPiece = type_of(pos.piece_on(bestMove.to_sq()));
        captureHistory[movedPiece][bestMove.to_sq()][capturedPiece] << captureBonus * S(1288, 1288) / S(1024, 1024);
    }

    if (prevSq != SQ_NONE && ((ss - 1)->moveCount == S(1, 1) + (ss - 1)->ttHit) && !pos.captured_piece())
        update_continuation_histories<PvNode>(ss - 1, pos.piece_on(prevSq), prevSq,
                                      -captureMalus * S(622, 622) / S(1024, 1024));

    for (Move move : capturesSearched)
    {
        movedPiece    = pos.moved_piece(move);
        capturedPiece = type_of(pos.piece_on(move.to_sq()));
        captureHistory[movedPiece][move.to_sq()][capturedPiece] << -captureMalus * S(1431, 1431) / S(1024, 1024);
    }
#undef S
}


// Updates histories of the move pairs
template<bool PvNode>
void update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus) {
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))
    static constexpr std::array<ConthistBonus, 6> conthist_bonuses = {
      {{S(1, 1), S(1108, 1108)}, {S(2, 2), S(652, 652)}, {S(3, 3), S(273, 273)}, {S(4, 4), S(572, 572)}, {S(5, 5), S(126, 126)}, {S(6, 6), S(449, 449)}}};

    for (const auto [i, weight] : conthist_bonuses)
    {
        if (ss->inCheck && i > S(2, 2)) break;
        if (((ss - i)->currentMove).is_ok())
            (*(ss - i)->continuationHistory)[pc][to] << (bonus * weight / S(1024, 1024)) + S(80, 80) * (i < S(2, 2));
    }
#undef S
}

// Updates move sorting heuristics
template<bool PvNode>
void update_quiet_histories(
  const Position& pos, Stack* ss, Search::Worker& workerThread, Move move, int bonus) {
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))
    Color us = pos.side_to_move();
    workerThread.mainHistory[us][move.from_to()] << bonus;

    if (ss->ply < LOW_PLY_HISTORY_SIZE)
        workerThread.lowPlyHistory[ss->ply][move.from_to()] << (bonus * S(771, 771) / S(1024, 1024)) + S(40, 40);

    update_continuation_histories<PvNode>(ss, pos.moved_piece(move), move.to_sq(),
                                  bonus * (bonus > S(0, 0) ? S(979, 979) : S(842, 842)) / S(1024, 1024));

    int pIndex = pawn_history_index(pos);
    workerThread.pawnHistory[pIndex][pos.moved_piece(move)][move.to_sq()]
      << (bonus * (bonus > S(0, 0) ? S(704, 704) : S(439, 439)) / S(1024, 1024)) + S(70, 70);
#undef S
}

}

// When playing with strength handicap
Move Skill::pick_best(const RootMoves& rootMoves, size_t multiPV) {
    static PRNG rng(now());
    Value  topScore = rootMoves[0].score;
    int    delta    = std::min(topScore - rootMoves[multiPV - 1].score, int(PawnValue));
    int    maxScore = -VALUE_INFINITE;
    double weakness = 120 - 2 * level;

    for (size_t i = 0; i < multiPV; ++i)
    {
        int push = int(weakness * int(topScore - rootMoves[i].score)
                       + delta * (rng.rand<unsigned>() % int(weakness))) / 128;
        if (rootMoves[i].score + push >= maxScore)
        {
            maxScore = rootMoves[i].score + push;
            best     = rootMoves[i].pv[0];
        }
    }
    return best;
}


// Used to print debug info and, more importantly, to detect
// when we are out of available time and thus stop the search.
void SearchManager::check_time(Search::Worker& worker) {
    if (--callsCnt > 0)
        return;

    callsCnt = worker.limits.nodes ? std::min(512, int(worker.limits.nodes / 1024)) : 512;

    static TimePoint lastInfoTime = now();

    TimePoint elapsed = tm.elapsed([&worker]() { return worker.threads.nodes_searched(); });
    TimePoint tick    = worker.limits.startTime + elapsed;

    if (tick - lastInfoTime >= 1000)
    {
        lastInfoTime = tick;
        dbg_print();
    }

    if (ponder)
        return;

    if (
      worker.completedDepth >= 1
      && ((worker.limits.use_time_management() && (elapsed > tm.maximum() || stopOnPonderhit))
          || (worker.limits.movetime && elapsed >= worker.limits.movetime)
          || (worker.limits.nodes && worker.threads.nodes_searched() >= worker.limits.nodes)))
        worker.threads.stop = worker.threads.abortedSearch = true;
}

// Used to correct and extend PVs for moves that have a TB score.
void syzygy_extend_pv(const OptionsMap&         options,
                      const Search::LimitsType& limits,
                      Position&                 pos,
                      RootMove&                 rootMove,
                      Value&                    v) {

    auto t_start      = std::chrono::steady_clock::now();
    int  moveOverhead = int(options["Move Overhead"]);
    bool rule50       = bool(options["Syzygy50MoveRule"]);

    auto time_abort = [&t_start, &moveOverhead, &limits]() -> bool {
        auto t_end = std::chrono::steady_clock::now();
        return limits.use_time_management()
            && 2 * std::chrono::duration<double, std::milli>(t_end - t_start).count()
                 > moveOverhead;
    };

    std::list<StateInfo> sts;
    auto& stRoot = sts.emplace_back();
    pos.do_move(rootMove.pv[0], stRoot);
    int ply = 1;

    while (size_t(ply) < rootMove.pv.size())
    {
        Move& pvMove = rootMove.pv[ply];
        RootMoves legalMoves;
        for (const auto& m : MoveList<LEGAL>(pos)) legalMoves.emplace_back(m);
        Tablebases::Config config = Tablebases::rank_root_moves(options, pos, legalMoves);
        RootMove&          rm     = *std::find(legalMoves.begin(), legalMoves.end(), pvMove);
        if (legalMoves[0].tbRank != rm.tbRank) break;
        ply++;
        auto& st = sts.emplace_back();
        pos.do_move(pvMove, st);
        if (config.rootInTB && ((rule50 && pos.is_draw(ply)) || pos.is_repetition(ply)))
        {
            pos.undo_move(pvMove);
            ply--;
            break;
        }
        if (config.rootInTB && time_abort()) break;
    }

    rootMove.pv.resize(ply);

    while (!(rule50 && pos.is_draw(0)))
    {
        if (time_abort()) break;
        RootMoves legalMoves;
        for (const auto& m : MoveList<LEGAL>(pos))
        {
            auto& rm = legalMoves.emplace_back(m);
            StateInfo tmpSI;
            pos.do_move(m, tmpSI);
            for (const auto& mOpp : MoveList<LEGAL>(pos))
                rm.tbRank -= pos.capture(mOpp) ? 100 : 1;
            pos.undo_move(m);
        }
        if (legalMoves.size() == 0) break;
        std::stable_sort(
          legalMoves.begin(), legalMoves.end(),
          [](const Search::RootMove& a, const Search::RootMove& b) { return a.tbRank > b.tbRank; });
        Tablebases::Config config = Tablebases::rank_root_moves(options, pos, legalMoves, true);
        if (!config.rootInTB || config.cardinality > 0) break;
        ply++;
        Move& pvMove = legalMoves[0].pv[0];
        rootMove.pv.push_back(pvMove);
        auto& st = sts.emplace_back();
        pos.do_move(pvMove, st);
    }

    if (pos.is_draw(0)) v = VALUE_DRAW;

    for (auto it = rootMove.pv.rbegin(); it != rootMove.pv.rend(); ++it)
        pos.undo_move(*it);

    if (time_abort())
        sync_cout << "info string Syzygy based PV extension requires more time, increase Move Overhead as needed." << sync_endl;
}

void SearchManager::pv(Search::Worker&           worker,
                       const ThreadPool&         threads,
                       const TranspositionTable& tt,
                       Depth                     depth) {

    const auto nodes     = threads.nodes_searched();
    auto&      rootMoves = worker.rootMoves;
    auto&      pos       = worker.rootPos;
    size_t     pvIdx     = worker.pvIdx;
    size_t     multiPV   = std::min(size_t(worker.options["MultiPV"]), rootMoves.size());
    uint64_t   tbHits    = threads.tb_hits() + (worker.tbConfig.rootInTB ? rootMoves.size() : 0);

    for (size_t i = 0; i < multiPV; ++i)
    {
        bool updated = rootMoves[i].score != -VALUE_INFINITE;
        if (depth == 1 && !updated && i > 0) continue;
        Depth d = updated ? depth : std::max(1, depth - 1);
        Value v = updated ? rootMoves[i].uciScore : rootMoves[i].previousScore;
        if (v == -VALUE_INFINITE) v = VALUE_ZERO;
        bool tb = worker.tbConfig.rootInTB && std::abs(v) <= VALUE_TB;
        v       = tb ? rootMoves[i].tbScore : v;
        bool isExact = i != pvIdx || tb || !updated;

        if (is_decisive(v) && std::abs(v) < VALUE_MATE_IN_MAX_PLY
            && ((!rootMoves[i].scoreLowerbound && !rootMoves[i].scoreUpperbound) || isExact))
            syzygy_extend_pv(worker.options, worker.limits, pos, rootMoves[i], v);

        std::string pv_str;
        for (Move m : rootMoves[i].pv) pv_str += UCIEngine::move(m, pos.is_chess960()) + " ";
        if (!pv_str.empty()) pv_str.pop_back();

        auto wdl   = worker.options["UCI_ShowWDL"] ? UCIEngine::wdl(v, pos) : "";
        auto bound = rootMoves[i].scoreLowerbound ? "lowerbound" : (rootMoves[i].scoreUpperbound ? "upperbound" : "");
        InfoFull info;
        info.depth    = d;
        info.selDepth = rootMoves[i].selDepth;
        info.multiPV  = i + 1;
        info.score    = {v, pos};
        info.wdl      = wdl;
        if (!isExact) info.bound = bound;
        TimePoint time = std::max(TimePoint(1), tm.elapsed_time());
        info.timeMs    = time;
        info.nodes     = nodes;
        info.nps       = nodes * 1000 / time;
        info.tbHits    = tbHits;
        info.pv        = pv_str;
        info.hashfull  = tt.hashfull();
        updates.onUpdateFull(info);
    }
}

bool RootMove::extract_ponder_from_tt(const TranspositionTable& tt, Position& pos) {
    StateInfo st;
    assert(pv.size() == 1);
    if (pv[0] == Move::none()) return false;
    pos.do_move(pv[0], st, &tt);
    auto [ttHit, ttData, ttWriter] = tt.probe(pos.key());
    if (ttHit)
    {
        if (MoveList<LEGAL>(pos).contains(ttData.move))
            pv.push_back(ttData.move);
    }
    pos.undo_move(pv[0]);
    return pv.size() > 1;
}

}  // namespace Stockfish
