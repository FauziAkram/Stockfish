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

template<bool PvNode> void update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus);
template<bool PvNode> void update_quiet_histories(const Position& pos, Stack* ss, Search::Worker& workerThread, Move move, int bonus);
template<bool PvNode> void update_all_stats(const Position& pos, Stack* ss, Search::Worker& workerThread, Move bestMove, Square prevSq, SearchedList& quietsSearched, SearchedList& capturesSearched, Depth depth, Move ttMove);


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
    workerThread.minorPieceCorrectionHistory[minor_piece_index(pos)][us] << bonus * S(153, 153) / 128;
    workerThread.nonPawnCorrectionHistory[non_pawn_index<WHITE>(pos)][WHITE][us]
      << bonus * S(nonPawnWeight, nonPawnWeight) / 128;
    workerThread.nonPawnCorrectionHistory[non_pawn_index<BLACK>(pos)][BLACK][us]
      << bonus * S(nonPawnWeight, nonPawnWeight) / 128;

    if (m.is_ok())
        (*(ss - 2)->continuationCorrectionHistory)[pos.piece_on(m.to_sq())][m.to_sq()]
          << bonus * S(153, 153) / 128;
#undef S
}

Value value_draw(size_t nodes) { return VALUE_DRAW - 1 + Value(nodes & 0x2); }
Value value_to_tt(Value v, int ply);
Value value_from_tt(Value v, int ply, int r50c);
void  update_pv(Move* pv, Move move, const Move* childPv);


}  // namespace

Search::Worker::Worker(SharedState&                    sharedState,
                       std::unique_ptr<ISearchManager> sm,
                       size_t                          threadId,
                       NumaReplicatedAccessToken       token) :
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
    (void) (networks[numaAccessToken]);
}

void Search::Worker::start_searching() {

    accumulatorStack.reset();

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
        threads.start_searching();
        iterative_deepening();
    }

    while (!threads.stop && (main_manager()->ponder || limits.infinite))
    {}

    threads.stop = true;
    threads.wait_for_search_finished();

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

    if (bestThread != this)
        main_manager()->pv(*bestThread, threads, tt, bestThread->completedDepth);

    std::string ponder;
    if (bestThread->rootMoves[0].pv.size() > 1
        || bestThread->rootMoves[0].extract_ponder_from_tt(tt, rootPos))
        ponder = UCIEngine::move(bestThread->rootMoves[0].pv[1], rootPos.is_chess960());

    auto bestmove = UCIEngine::move(bestThread->rootMoves[0].pv[0], rootPos.is_chess960());
    main_manager()->updates.onBestmove(bestmove, ponder);
}

void Search::Worker::iterative_deepening() {

    SearchManager* mainThread = (is_mainthread() ? main_manager() : nullptr);

    Move pv[MAX_PLY + 1];

    Depth lastBestMoveDepth = 0;
    Value lastBestScore     = -VALUE_INFINITE;
    auto  lastBestPV        = std::vector{Move::none()};

    Value  alpha, beta;
    Value  bestValue     = -VALUE_INFINITE;
    Color  us            = rootPos.side_to_move();
    double timeReduction = 1, totBestMoveChanges = 0;
    int    delta, iterIdx                        = 0;

    Stack  stack[MAX_PLY + 10] = {};
    Stack* ss                  = stack + 7;

    for (int i = 7; i > 0; --i)
    {
        (ss - i)->continuationHistory =
          &continuationHistory[0][0][NO_PIECE][0];
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

    if (skill.enabled())
        multiPV = std::max(multiPV, size_t(4));

    multiPV = std::min(multiPV, rootMoves.size());
    int searchAgainCounter = 0;
    lowPlyHistory.fill(89);

    while (++rootDepth < MAX_PLY && !threads.stop
           && !(limits.depth && mainThread && rootDepth > limits.depth))
    {
        if (mainThread)
            totBestMoveChanges /= 2;

        for (RootMove& rm : rootMoves)
            rm.previousScore = rm.score;

        size_t pvFirst = 0;
        pvLast         = 0;

        if (!threads.increaseDepth)
            searchAgainCounter++;

        for (pvIdx = 0; pvIdx < multiPV; ++pvIdx)
        {
            if (pvIdx == pvLast)
            {
                pvFirst = pvLast;
                for (pvLast++; pvLast < rootMoves.size(); pvLast++)
                    if (rootMoves[pvLast].tbRank != rootMoves[pvFirst].tbRank)
                        break;
            }

            selDepth = 0;
            delta     = 5 + std::abs(rootMoves[pvIdx].meanSquaredScore) / 11131;
            Value avg = rootMoves[pvIdx].averageScore;
            alpha     = std::max(avg - delta, -VALUE_INFINITE);
            beta      = std::min(avg + delta, VALUE_INFINITE);
            optimism[us]  = 136 * avg / (std::abs(avg) + 93);
            optimism[~us] = -optimism[us];

            int failedHighCnt = 0;
            while (true)
            {
                Depth adjustedDepth =
                  std::max(1, rootDepth - failedHighCnt - 3 * (searchAgainCounter + 1) / 4);
                rootDelta = beta - alpha;
                bestValue = search<Root>(rootPos, ss, alpha, beta, adjustedDepth, false);

                std::stable_sort(rootMoves.begin() + pvIdx, rootMoves.begin() + pvLast);

                if (threads.stop)
                    break;

                if (mainThread && multiPV == 1 && (bestValue <= alpha || bestValue >= beta)
                    && nodes > 10000000)
                    main_manager()->pv(*this, threads, tt, rootDepth);

                if (bestValue <= alpha)
                {
                    beta  = (3 * alpha + beta) / 4;
                    alpha = std::max(bestValue - delta, -VALUE_INFINITE);
                    failedHighCnt = 0;
                    if (mainThread) mainThread->stopOnPonderhit = false;
                }
                else if (bestValue >= beta)
                {
                    beta = std::min(bestValue + delta, VALUE_INFINITE);
                    ++failedHighCnt;
                }
                else
                    break;

                delta += delta / 3;

                assert(alpha >= -VALUE_INFINITE && beta <= VALUE_INFINITE);
            }

            std::stable_sort(rootMoves.begin() + pvFirst, rootMoves.begin() + pvIdx + 1);

            if (mainThread
                && (threads.stop || pvIdx + 1 == multiPV || nodes > 10000000)
                && !(threads.abortedSearch && is_loss(rootMoves[0].uciScore)))
                main_manager()->pv(*this, threads, tt, rootDepth);

            if (threads.stop)
                break;
        }

        if (!threads.stop)
            completedDepth = rootDepth;

        if (threads.abortedSearch && rootMoves[0].score != -VALUE_INFINITE
            && is_loss(rootMoves[0].score))
        {
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

        if (limits.mate && rootMoves[0].score == rootMoves[0].uciScore
            && ((rootMoves[0].score >= VALUE_MATE_IN_MAX_PLY
                 && VALUE_MATE - rootMoves[0].score <= 2 * limits.mate)
                || (rootMoves[0].score != -VALUE_INFINITE
                    && rootMoves[0].score <= VALUE_MATED_IN_MAX_PLY
                    && VALUE_MATE + rootMoves[0].score <= 2 * limits.mate)))
            threads.stop = true;

        if (skill.enabled() && skill.time_to_pick(rootDepth))
            skill.pick_best(rootMoves, multiPV);

        for (auto&& th : threads)
        {
            totBestMoveChanges += th->worker->bestMoveChanges;
            th->worker->bestMoveChanges = 0;
        }

        if (limits.use_time_management() && !threads.stop && !mainThread->stopOnPonderhit)
        {
            uint64_t nodesEffort =
              rootMoves[0].effort * 100000 / std::max(size_t(1), size_t(nodes));
            double fallingEval =
              (11.396 + 2.035 * (mainThread->bestPreviousAverageScore - bestValue)
               + 0.968 * (mainThread->iterValue[iterIdx] - bestValue)) / 100.0;
            fallingEval = std::clamp(fallingEval, 0.5786, 1.6752);
            double k      = 0.527;
            double center = lastBestMoveDepth + 11.0;
            timeReduction = 0.8 + 0.84 / (1.077 + std::exp(-k * (completedDepth - center)));
            double reduction =
              (1.4540 + mainThread->previousTimeReduction) / (2.1593 * timeReduction);
            double bestMoveInstability = 0.9929 + 1.8519 * totBestMoveChanges / threads.size();
            double totalTime = mainThread->tm.optimum() * fallingEval * reduction * bestMoveInstability;
            if (rootMoves.size() == 1) totalTime = std::min(500.0, totalTime);
            auto elapsedTime = elapsed();
            if (completedDepth >= 10 && nodesEffort >= 97056 && elapsedTime > totalTime * 0.6540 && !mainThread->ponder)
                threads.stop = true;
            if (elapsedTime > std::min(totalTime, double(mainThread->tm.maximum())))
            {
                if (mainThread->ponder) mainThread->stopOnPonderhit = true;
                else threads.stop = true;
            }
            else
                threads.increaseDepth = mainThread->ponder || elapsedTime <= totalTime * 0.5138;
        }

        mainThread->iterValue[iterIdx] = bestValue;
        iterIdx                        = (iterIdx + 1) & 3;
    }

    if (!mainThread)
        return;

    mainThread->previousTimeReduction = timeReduction;

    if (skill.enabled())
        std::swap(rootMoves[0],
                  *std::find(rootMoves.begin(), rootMoves.end(),
                             skill.best ? skill.best : skill.pick_best(rootMoves, multiPV)));
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
void Search::Worker::undo_move(Position& pos, const Move move) { pos.undo_move(move); accumulatorStack.pop(); }
void Search::Worker::undo_null_move(Position& pos) { pos.undo_null_move(); }


void Search::Worker::clear() {
    mainHistory.fill(64);
    captureHistory.fill(-753);
    pawnHistory.fill(-1275);
    pawnCorrectionHistory.fill(5);
    minorPieceCorrectionHistory.fill(0);
    nonPawnCorrectionHistory.fill(0);
    ttMoveHistory = 0;
    for (auto& to : continuationCorrectionHistory) for (auto& h : to) h.fill(8);
    for (bool inCheck : {false, true})
        for (StatsType c : {NoCaptures, Captures})
            for (auto& to : continuationHistory[inCheck][c])
                for (auto& h : to) h.fill(-494);
    for (size_t i = 1; i < reductions.size(); ++i)
        reductions[i] = int(2782 / 128.0 * std::log(i));
    refreshTable.clear(networks[numaAccessToken]);
}


template<NodeType nodeType>
Value Search::Worker::search(
  Position& pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode) {

    constexpr bool PvNode   = nodeType != NonPV;
    constexpr bool rootNode = nodeType == Root;
    const bool     allNode  = !(PvNode || cutNode);

#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))

    if (depth <= 0)
    {
        constexpr auto nt = PvNode ? PV : NonPV;
        return qsearch<nt>(pos, ss, alpha, beta);
    }

    depth = std::min(depth, Depth(MAX_PLY - 1));
    if (!rootNode && alpha < VALUE_DRAW && pos.upcoming_repetition(ss->ply))
    {
        alpha = value_draw(nodes);
        if (alpha >= beta) return alpha;
    }

    assert(-VALUE_INFINITE <= alpha && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - 1));
    assert(0 < depth && depth < MAX_PLY);
    assert(!(PvNode && cutNode));

    Move pv[MAX_PLY + 1];
    StateInfo st;
    Key posKey;
    Move move, excludedMove, bestMove;
    Depth extension, newDepth;
    Value bestValue, value, eval, maxValue, probCutBeta;
    bool  givesCheck, improving, priorCapture, opponentWorsening;
    bool  capture, ttCapture;
    int   priorReduction;
    Piece movedPiece;
    SearchedList capturesSearched, quietsSearched;

    ss->inCheck   = pos.checkers();
    priorCapture  = pos.captured_piece();
    Color us      = pos.side_to_move();
    ss->moveCount = 0;
    bestValue     = -VALUE_INFINITE;
    maxValue      = VALUE_INFINITE;

    if (is_mainthread()) main_manager()->check_time(*this);
    if (PvNode && selDepth < ss->ply + 1) selDepth = ss->ply + 1;

    if (!rootNode)
    {
        if (threads.stop.load(std::memory_order_relaxed) || pos.is_draw(ss->ply) || ss->ply >= MAX_PLY)
            return (ss->ply >= MAX_PLY && !ss->inCheck) ? evaluate(pos) : value_draw(nodes);
        alpha = std::max(mated_in(ss->ply), alpha);
        beta  = std::min(mate_in(ss->ply + 1), beta);
        if (alpha >= beta) return alpha;
    }

    assert(0 <= ss->ply && ss->ply < MAX_PLY);

    Square prevSq  = ((ss - 1)->currentMove).is_ok() ? ((ss - 1)->currentMove).to_sq() : SQ_NONE;
    bestMove       = Move::none();
    priorReduction = (ss - 1)->reduction;
    (ss - 1)->reduction = 0;
    ss->statScore       = 0;
    (ss + 2)->cutoffCnt = 0;

    excludedMove                   = ss->excludedMove;
    posKey                         = pos.key();
    auto [ttHit, ttData, ttWriter] = tt.probe(posKey);
    ss->ttHit    = ttHit;
    ttData.move  = rootNode ? rootMoves[pvIdx].pv[0] : ttHit ? ttData.move : Move::none();
    ttData.value = ttHit ? value_from_tt(ttData.value, ss->ply, pos.rule50_count()) : VALUE_NONE;
    ss->ttPv     = excludedMove ? ss->ttPv : PvNode || (ttHit && ttData.is_pv);
    ttCapture    = ttData.move && pos.capture_stage(ttData.move);

    if (!PvNode && !excludedMove && ttData.depth > depth - 1 * (ttData.value <= beta)
        && is_valid(ttData.value)
        && (ttData.bound & (ttData.value >= beta ? BOUND_LOWER : BOUND_UPPER))
        && (cutNode == (ttData.value >= beta) || depth > 5))
    {
        if (ttData.move && ttData.value >= beta)
        {
            if (!ttCapture)
                update_quiet_histories<PvNode>(pos, ss, *this, ttData.move,
                                       std::min(S(127, 127) * depth - S(74, 74), S(1063, 1063)));
            if (prevSq != SQ_NONE && (ss - 1)->moveCount <= 3 && !priorCapture)
                update_continuation_histories<PvNode>(ss - 1, pos.piece_on(prevSq), prevSq, S(-2128, -2128));
        }
        if (pos.rule50_count() < 91)
        {
            if (depth >= 8 && ttData.move && pos.pseudo_legal(ttData.move) && pos.legal(ttData.move)
                && !is_decisive(ttData.value))
            {
                pos.do_move(ttData.move, st);
                auto [ttHitNext, ttDataNext, ttWriterNext] = tt.probe(pos.key());
                pos.undo_move(ttData.move);
                if (!is_valid(ttDataNext.value)) return ttData.value;
                if ((ttData.value >= beta) == (-ttDataNext.value >= beta)) return ttData.value;
            }
            else return ttData.value;
        }
    }

    if (!rootNode && !excludedMove && tbConfig.cardinality)
    {
        int piecesCount = pos.count<ALL_PIECES>();
        if (piecesCount <= tbConfig.cardinality && (piecesCount < tbConfig.cardinality || depth >= tbConfig.probeDepth)
            && pos.rule50_count() == 0 && !pos.can_castle(ANY_CASTLING))
        {
            TB::ProbeState err;
            TB::WDLScore wdl = Tablebases::probe_wdl(pos, &err);
            if (is_mainthread()) main_manager()->callsCnt = 0;
            if (err != TB::ProbeState::FAIL)
            {
                tbHits.fetch_add(1, std::memory_order_relaxed);
                int drawScore = tbConfig.useRule50 ? 1 : 0;
                Value tbValue = VALUE_TB - ss->ply;
                value = wdl < -drawScore ? -tbValue : wdl > drawScore  ? tbValue : VALUE_DRAW + 2 * wdl * drawScore;
                Bound b = wdl < -drawScore ? BOUND_UPPER : wdl > drawScore  ? BOUND_LOWER : BOUND_EXACT;
                if (b == BOUND_EXACT || (b == BOUND_LOWER ? value >= beta : value <= alpha))
                {
                    ttWriter.write(posKey, value_to_tt(value, ss->ply), ss->ttPv, b,
                                   std::min(Depth(MAX_PLY - 1), depth + 6), Move::none(), VALUE_NONE, tt.generation());
                    return value;
                }
                if (PvNode)
                {
                    if (b == BOUND_LOWER) bestValue = value, alpha = std::max(alpha, bestValue);
                    else maxValue = value;
                }
            }
        }
    }

    Value unadjustedStaticEval = VALUE_NONE;
    const auto correctionValue = correction_value<PvNode>(*this, pos, ss);
    if (ss->inCheck)
    {
        ss->staticEval = eval = (ss - 2)->staticEval;
        improving             = false;
        goto moves_loop;
    }
    else if (excludedMove) unadjustedStaticEval = eval = ss->staticEval;
    else if (ss->ttHit)
    {
        unadjustedStaticEval = ttData.eval;
        if (!is_valid(unadjustedStaticEval)) unadjustedStaticEval = evaluate(pos);
        ss->staticEval = eval = to_corrected_static_eval(unadjustedStaticEval, correctionValue);
        if (is_valid(ttData.value) && (ttData.bound & (ttData.value > eval ? BOUND_LOWER : BOUND_UPPER)))
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
        int bonus = std::clamp(-10 * int((ss - 1)->staticEval + ss->staticEval), -1979, 1561) + 630;
        mainHistory[~us][((ss - 1)->currentMove).from_to()] << bonus * 935 / 1024;
        if (!ttHit && type_of(pos.piece_on(prevSq)) != PAWN && ((ss - 1)->currentMove).type_of() != PROMOTION)
            pawnHistory[pawn_history_index(pos)][pos.piece_on(prevSq)][prevSq] << bonus * 1428 / 1024;
    }

    improving         = ss->staticEval > (ss - 2)->staticEval;
    opponentWorsening = ss->staticEval > -(ss - 1)->staticEval;

    if (priorReduction >= (depth < 10 ? 1 : 3) && !opponentWorsening) depth++;
    if (priorReduction >= 2 && depth >= 2 && ss->staticEval + (ss - 1)->staticEval > 177) depth--;

    if (!PvNode && eval < alpha - 495 - 290 * depth * depth)
        return qsearch<NonPV>(pos, ss, alpha, beta);

    {
        auto futility_margin = [&](Depth d) {
            Value futilityMult = 90 - 20 * (cutNode && !ss->ttHit);
            return futilityMult * d - improving * futilityMult * 2 - opponentWorsening * futilityMult / 3
                 + (ss - 1)->statScore / 356 + std::abs(correctionValue) / 171290;
        };
        if (!ss->ttPv && depth < 14 && eval - futility_margin(depth) >= beta && eval >= beta
            && (!ttData.move || ttCapture) && !is_loss(beta) && !is_win(eval))
            return beta + (eval - beta) / 3;
    }

    if (cutNode && ss->staticEval >= beta - 19 * depth + 403 && !excludedMove
        && pos.non_pawn_material(us) && ss->ply >= nmpMinPly && !is_loss(beta))
    {
        assert((ss - 1)->currentMove != Move::null());
        Depth R = 7 + depth / 3;
        ss->currentMove = Move::null();
        ss->continuationHistory = &continuationHistory[0][0][NO_PIECE][0];
        ss->continuationCorrectionHistory = &continuationCorrectionHistory[NO_PIECE][0];
        do_null_move(pos, st);
        Value nullValue = -search<NonPV>(pos, ss + 1, -beta, -beta + 1, depth - R, false);
        undo_null_move(pos);
        if (nullValue >= beta && !is_win(nullValue))
        {
            if (nmpMinPly || depth < 16) return nullValue;
            assert(!nmpMinPly);
            nmpMinPly = ss->ply + 3 * (depth - R) / 4;
            Value v = search<NonPV>(pos, ss, beta - 1, beta, depth - R, false);
            nmpMinPly = 0;
            if (v >= beta) return nullValue;
        }
    }

    improving |= ss->staticEval >= beta;
    if (!allNode && depth >= 6 && !ttData.move && priorReduction <= 3) depth--;

    probCutBeta = beta + 215 - 60 * improving;
    if (depth >= 3 && !is_decisive(beta) && !(is_valid(ttData.value) && ttData.value < probCutBeta))
    {
        assert(probCutBeta < VALUE_INFINITE && probCutBeta > beta);
        MovePicker mp(pos, ttData.move, probCutBeta - ss->staticEval, &captureHistory);
        Depth dynamicReduction = (ss->staticEval - beta) / 300;
        Depth probCutDepth     = std::max(Depth(0), depth - 5 - dynamicReduction);
        while ((move = mp.next_move()) != Move::none())
        {
            assert(move.is_ok());
            if (move == excludedMove || !pos.legal(move)) continue;
            assert(pos.capture_stage(move));
            movedPiece = pos.moved_piece(move);
            do_move(pos, move, st, ss);
            value = -qsearch<NonPV>(pos, ss + 1, -probCutBeta, -probCutBeta + 1);
            if (value >= probCutBeta && probCutDepth > 0)
                value = -search<NonPV>(pos, ss + 1, -probCutBeta, -probCutBeta + 1, probCutDepth, !cutNode);
            undo_move(pos, move);
            if (value >= probCutBeta)
            {
                ttWriter.write(posKey, value_to_tt(value, ss->ply), ss->ttPv, BOUND_LOWER,
                               probCutDepth + 1, move, unadjustedStaticEval, tt.generation());
                if (!is_decisive(value)) return value - (probCutBeta - beta);
            }
        }
    }

moves_loop:
    probCutBeta = beta + 417;
    if ((ttData.bound & BOUND_LOWER) && ttData.depth >= depth - 4 && ttData.value >= probCutBeta
        && !is_decisive(beta) && is_valid(ttData.value) && !is_decisive(ttData.value))
        return probCutBeta;

    const PieceToHistory* contHist[] = {
      (ss - 1)->continuationHistory, (ss - 2)->continuationHistory, (ss - 3)->continuationHistory,
      (ss - 4)->continuationHistory, (ss - 5)->continuationHistory, (ss - 6)->continuationHistory};
    MovePicker mp(pos, ttData.move, depth, &mainHistory, &lowPlyHistory, &captureHistory, contHist, &pawnHistory, ss->ply);

    value = bestValue;
    int moveCount = 0;

    while ((move = mp.next_move()) != Move::none())
    {
        assert(move.is_ok());
        if (move == excludedMove) continue;
        if (!pos.legal(move)) continue;
        if (rootNode && !std::count(rootMoves.begin() + pvIdx, rootMoves.begin() + pvLast, move)) continue;

        ss->moveCount = ++moveCount;
        if (rootNode && is_mainthread() && nodes > 10000000)
            main_manager()->updates.onIter({depth, UCIEngine::move(move, pos.is_chess960()), moveCount + pvIdx});
        if (PvNode) (ss + 1)->pv = nullptr;

        extension  = 0;
        capture    = pos.capture_stage(move);
        movedPiece = pos.moved_piece(move);
        givesCheck = pos.gives_check(move);
        (ss + 1)->quietMoveStreak = (!capture && !givesCheck) ? (ss->quietMoveStreak + 1) : 0;
        newDepth = depth - 1;
        int delta = beta - alpha;
        Depth r = this->reduction<PvNode>(improving, depth, moveCount, delta);
        if (ss->ttPv) r += 931;

        if (!rootNode && pos.non_pawn_material(us) && !is_loss(bestValue))
        {
            if (moveCount >= (3 + depth * depth) / (2 - improving)) mp.skip_quiet_moves();
            int lmrDepth = newDepth - r / 1024;
            if (capture || givesCheck)
            {
                Piece capturedPiece = pos.piece_on(move.to_sq());
                int captHist = captureHistory[movedPiece][move.to_sq()][type_of(capturedPiece)];
                if (!givesCheck && lmrDepth < 7 && !ss->inCheck)
                {
                    Value futilityValue = ss->staticEval + 232 + 224 * lmrDepth + PieceValue[capturedPiece] + 131 * captHist / 1024;
                    if (futilityValue <= alpha) continue;
                }
                int margin = std::clamp(158 * depth + captHist / 31, 0, 283 * depth);
                if (!pos.see_ge(move, -margin))
                {
                    bool mayStalemateTrap = depth > 2 && alpha < 0 && pos.non_pawn_material(us) == PieceValue[movedPiece]
                                         && PieceValue[movedPiece] >= RookValue
                                         && !(attacks_bb<KING>(pos.square<KING>(us)) & move.from_sq())
                                         && !mp.can_move_king_or_pawn();
                    if (!mayStalemateTrap) continue;
                }
            }
            else
            {
                int history = (*contHist[0])[movedPiece][move.to_sq()] + (*contHist[1])[movedPiece][move.to_sq()]
                            + pawnHistory[pawn_history_index(pos)][movedPiece][move.to_sq()];
                if (history < -4361 * depth) continue;
                history += 71 * mainHistory[us][move.from_to()] / 32;
                lmrDepth += history / 3233;
                Value baseFutility = (bestMove ? 46 : 230);
                Value futilityValue = ss->staticEval + baseFutility + 131 * lmrDepth + 91 * (ss->staticEval > alpha);
                if (!ss->inCheck && lmrDepth < 11 && futilityValue <= alpha)
                {
                    if (bestValue <= futilityValue && !is_decisive(bestValue) && !is_win(futilityValue))
                        bestValue = futilityValue;
                    continue;
                }
                lmrDepth = std::max(lmrDepth, 0);
                if (!pos.see_ge(move, -26 * lmrDepth * lmrDepth)) continue;
            }
        }

        if (!rootNode && move == ttData.move && !excludedMove
            && depth >= S(6, 6) - (completedDepth > 26) + ss->ttPv && is_valid(ttData.value)
            && !is_decisive(ttData.value) && (ttData.bound & BOUND_LOWER)
            && ttData.depth >= depth - 3)
        {
            Value singularBeta  = ttData.value - (S(56, 56) + S(79, 79) * (ss->ttPv && !PvNode)) * depth / 58;
            Depth singularDepth = newDepth / 2;
            ss->excludedMove = move;
            value = search<NonPV>(pos, ss, singularBeta - 1, singularBeta, singularDepth, cutNode);
            ss->excludedMove = Move::none();
            if (value < singularBeta)
            {
                int corrValAdj   = std::abs(correctionValue) / S(249096, 249096);
                int doubleMargin = S(4, 209) - S(223, 223) * !ttCapture - corrValAdj
                                 - S(959, 959) * ttMoveHistory / 131072 - (ss->ply > rootDepth) * S(45, 45);
                int tripleMargin = S(80, 356) - S(249, 249) * !ttCapture + S(86, 86) * ss->ttPv - corrValAdj
                                 - (ss->ply * 2 > rootDepth * 3) * S(53, 53);
                extension = 1 + (value < singularBeta - doubleMargin) + (value < singularBeta - tripleMargin);
                depth++;
            }
            else if (value >= beta && !is_decisive(value)) return value;
            else if (ttData.value >= beta) extension = -3;
            else if (cutNode) extension = -2;
        }

        do_move(pos, move, st, givesCheck, ss);
        newDepth += extension;
        uint64_t nodeCount = rootNode ? uint64_t(nodes) : 0;

        if (ss->ttPv)
            r -= S(2510, 3473) + (ttData.value > alpha) * S(916, 916)
               + (ttData.depth >= depth) * (S(943, 943) + cutNode * S(1180, 1180));
        
        r += 650;
        r -= moveCount * 69;
        r -= std::abs(correctionValue) / 27160;
        if (cutNode) r += 3000 + 1024 * !ttData.move;
        if (ttCapture) r += 1350;
        if ((ss + 1)->cutoffCnt > 2) r += 935 + allNode * 763;
        r += (ss + 1)->quietMoveStreak * 51;
        if (move == ttData.move) r -= 2043;

        if (capture)
            ss->statScore = 782 * int(PieceValue[pos.captured_piece()]) / 128
                          + captureHistory[movedPiece][move.to_sq()][type_of(pos.captured_piece())];
        else
            ss->statScore = 2 * mainHistory[us][move.from_to()]
                          + (*contHist[0])[movedPiece][move.to_sq()]
                          + (*contHist[1])[movedPiece][move.to_sq()];
        r -= ss->statScore * 789 / 8192;

        if (depth >= 2 && moveCount > 1)
        {
            Depth d = std::max(Depth(1), std::min(newDepth - r / 1024, newDepth + S(1, 2))) + S(0, 1);
            ss->reduction = newDepth - d;
            value         = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, d, true);
            ss->reduction = 0;
            if (value > alpha)
            {
                const bool doDeeperSearch = d < newDepth && value > (bestValue + 43 + 2 * newDepth);
                const bool doShallowerSearch = value < bestValue + 9;
                newDepth += doDeeperSearch - doShallowerSearch;
                if (newDepth > d)
                    value = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode);
                update_continuation_histories<PvNode>(ss, movedPiece, move.to_sq(), S(1412, 1412));
            }
        }
        else if (!PvNode || moveCount > 1)
        {
            if (!ttData.move) r += 1139;
            const int threshold1 = depth <= 4 ? 2000 : 3200;
            const int threshold2 = depth <= 4 ? 3500 : 4600;
            value = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha,
                                   newDepth - (r > threshold1) - (r > threshold2 && newDepth > 2), !cutNode);
        }

        if (PvNode && (moveCount == 1 || value > alpha))
        {
            (ss + 1)->pv    = pv;
            (ss + 1)->pv[0] = Move::none();
            if (move == ttData.move && rootDepth > 8) newDepth = std::max(newDepth, Depth(1));
            value = -search<PV>(pos, ss + 1, -beta, -alpha, newDepth, false);
        }

        undo_move(pos, move);
        assert(value > -VALUE_INFINITE && value < VALUE_INFINITE);

        if (threads.stop.load(std::memory_order_relaxed)) return VALUE_ZERO;

        if (rootNode)
        {
            RootMove& rm = *std::find(rootMoves.begin(), rootMoves.end(), move);
            rm.effort += nodes - nodeCount;
            rm.averageScore = rm.averageScore != -VALUE_INFINITE ? (value + rm.averageScore) / 2 : value;
            rm.meanSquaredScore = rm.meanSquaredScore != -VALUE_INFINITE * VALUE_INFINITE
                                  ? (value * std::abs(value) + rm.meanSquaredScore) / 2 : value * std::abs(value);
            if (moveCount == 1 || value > alpha)
            {
                rm.score = rm.uciScore = value;
                rm.selDepth            = selDepth;
                rm.scoreLowerbound = rm.scoreUpperbound = false;
                if (value >= beta) { rm.scoreLowerbound = true; rm.uciScore = beta; }
                else if (value <= alpha) { rm.scoreUpperbound = true; rm.uciScore = alpha; }
                rm.pv.resize(1);
                assert((ss + 1)->pv);
                for (Move* m = (ss + 1)->pv; *m != Move::none(); ++m) rm.pv.push_back(*m);
                if (moveCount > 1 && !pvIdx) ++bestMoveChanges;
            } else rm.score = -VALUE_INFINITE;
        }

        int inc = (value == bestValue && ss->ply + 2 >= rootDepth && (int(nodes) & 14) == 0 && !is_win(std::abs(value) + 1));
        if (value + inc > bestValue)
        {
            bestValue = value;
            if (value + inc > alpha)
            {
                bestMove = move;
                if (PvNode && !rootNode) update_pv(ss->pv, move, (ss + 1)->pv);
                if (value >= beta)
                {
                    ss->cutoffCnt += (extension < 2) || PvNode;
                    assert(value >= beta);
                    break;
                }
                if (depth > 2 && depth < 16 && !is_decisive(value)) depth -= 2;
                assert(depth > 0);
                alpha = value;
            }
        }
        if (move != bestMove && moveCount <= SEARCHEDLIST_CAPACITY)
        {
            if (capture) capturesSearched.push_back(move);
            else quietsSearched.push_back(move);
        }
    }

    assert(moveCount || !ss->inCheck || excludedMove || !MoveList<LEGAL>(pos).size());
    if (bestValue >= beta && !is_decisive(bestValue) && !is_decisive(alpha))
        bestValue = (bestValue * depth + beta) / (depth + 1);

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
        int bonusScale = -215;
        bonusScale += std::min(-(ss - 1)->statScore / 103, 337);
        bonusScale += std::min(64 * depth, 552);
        bonusScale += 177 * ((ss - 1)->moveCount > 8);
        bonusScale += 141 * (!ss->inCheck && bestValue <= ss->staticEval - 94);
        bonusScale += 141 * (!(ss - 1)->inCheck && bestValue <= -(ss - 1)->staticEval - 76);
        bonusScale = std::max(bonusScale, 0);
        const int scaledBonus = std::min(155 * depth - 88, 1416) * bonusScale;
        update_continuation_histories<PvNode>(ss - 1, pos.piece_on(prevSq), prevSq, scaledBonus * S(397, 397) / 32768);
        mainHistory[~us][((ss - 1)->currentMove).from_to()] << scaledBonus * 224 / 32768;
        if (type_of(pos.piece_on(prevSq)) != PAWN && ((ss - 1)->currentMove).type_of() != PROMOTION)
            pawnHistory[pawn_history_index(pos)][pos.piece_on(prevSq)][prevSq] << scaledBonus * 1127 / 32768;
    }
    else if (priorCapture && prevSq != SQ_NONE)
    {
        Piece capturedPiece = pos.captured_piece();
        assert(capturedPiece != NO_PIECE);
        captureHistory[pos.piece_on(prevSq)][prevSq][type_of(capturedPiece)] << 1042;
    }

    if (PvNode) bestValue = std::min(bestValue, maxValue);
    if (bestValue <= alpha) ss->ttPv = ss->ttPv || (ss - 1)->ttPv;

    if (!excludedMove && !(rootNode && pvIdx))
        ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), ss->ttPv,
                       bestValue >= beta    ? BOUND_LOWER
                       : PvNode && bestMove ? BOUND_EXACT
                                            : BOUND_UPPER,
                       moveCount != 0 ? depth : std::min(Depth(MAX_PLY - 1), depth + 6), bestMove,
                       unadjustedStaticEval, tt.generation());

    if (!ss->inCheck && !(bestMove && pos.capture(bestMove))
        && ((bestValue < ss->staticEval && bestValue < beta)
            || (bestValue > ss->staticEval && bestMove)))
    {
        auto bonus = std::clamp(int(bestValue - ss->staticEval) * depth / 8, -CORRECTION_HISTORY_LIMIT / 4, CORRECTION_HISTORY_LIMIT / 4);
        update_correction_history<PvNode>(pos, ss, *this, bonus);
    }

    assert(bestValue > -VALUE_INFINITE && bestValue < VALUE_INFINITE);
#undef S
    return bestValue;
}


template<NodeType nodeType>
Value Search::Worker::qsearch(Position& pos, Stack* ss, Value alpha, Value beta) {

    static_assert(nodeType != Root);
    constexpr bool PvNode = nodeType == PV;
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))

    assert(alpha >= -VALUE_INFINITE && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - 1));

    if (alpha < VALUE_DRAW && pos.upcoming_repetition(ss->ply))
    {
        alpha = value_draw(nodes);
        if (alpha >= beta) return alpha;
    }

    Move pv[MAX_PLY + 1];
    StateInfo st;
    Key posKey;
    Move move, bestMove;
    Value bestValue, value, futilityBase;
    bool pvHit, givesCheck, capture;
    int moveCount;

    if (PvNode) { (ss + 1)->pv = pv; ss->pv[0] = Move::none(); }
    bestMove    = Move::none();
    ss->inCheck = pos.checkers();
    moveCount   = 0;

    if (PvNode && selDepth < ss->ply + 1) selDepth = ss->ply + 1;

    if (pos.is_draw(ss->ply) || ss->ply >= MAX_PLY)
        return (ss->ply >= MAX_PLY && !ss->inCheck) ? evaluate(pos) : VALUE_DRAW;

    assert(0 <= ss->ply && ss->ply < MAX_PLY);

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
    if (ss->inCheck) bestValue = futilityBase = -VALUE_INFINITE;
    else
    {
        const auto correctionValue = correction_value<PvNode>(*this, pos, ss);
        if (ss->ttHit)
        {
            unadjustedStaticEval = ttData.eval;
            if (!is_valid(unadjustedStaticEval)) unadjustedStaticEval = evaluate(pos);
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
            if (!is_decisive(bestValue)) bestValue = (bestValue + beta) / 2;
            if (!ss->ttHit)
                ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), false, BOUND_LOWER,
                               DEPTH_UNSEARCHED, Move::none(), unadjustedStaticEval, tt.generation());
            return bestValue;
        }
        if (bestValue > alpha) alpha = bestValue;
        futilityBase = ss->staticEval + S(352, 352);
    }

    const PieceToHistory* contHist[] = {(ss - 1)->continuationHistory, (ss - 2)->continuationHistory};
    Square prevSq = ((ss - 1)->currentMove).is_ok() ? ((ss - 1)->currentMove).to_sq() : SQ_NONE;

    MovePicker mp(pos, ttData.move, DEPTH_QS, &mainHistory, &lowPlyHistory, &captureHistory, contHist, &pawnHistory, ss->ply);

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
                if (moveCount > 2) continue;
                Value futilityValue = futilityBase + PieceValue[pos.piece_on(move.to_sq())];
                if (futilityValue <= alpha) { bestValue = std::max(bestValue, futilityValue); continue; }
                if (!pos.see_ge(move, alpha - futilityBase)) { bestValue = std::min(alpha, futilityBase); continue; }
            }
            if (!capture && (*contHist[0])[pos.moved_piece(move)][move.to_sq()]
                   + pawnHistory[pawn_history_index(pos)][pos.moved_piece(move)][move.to_sq()] <= S(5868, 5868))
                continue;
            if (!pos.see_ge(move, -S(74, 74))) continue;
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
    if (!is_decisive(bestValue) && bestValue > beta) bestValue = (bestValue + beta) / 2;

    Color us = pos.side_to_move();
    if (!ss->inCheck && !moveCount && !pos.non_pawn_material(us) && type_of(pos.captured_piece()) >= ROOK)
    {
        if (!((us == WHITE ? shift<NORTH>(pos.pieces(us, PAWN)) : shift<SOUTH>(pos.pieces(us, PAWN))) & ~pos.pieces()))
        {
            pos.state()->checkersBB = Rank1BB;
            if (!MoveList<LEGAL>(pos).size()) bestValue = VALUE_DRAW;
            pos.state()->checkersBB = 0;
        }
    }

    ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), pvHit,
                   bestValue >= beta ? BOUND_LOWER : BOUND_UPPER, DEPTH_QS, bestMove,
                   unadjustedStaticEval, tt.generation());

    assert(bestValue > -VALUE_INFINITE && bestValue < VALUE_INFINITE);
#undef S
    return bestValue;
}

template<bool PvNode>
Depth Search::Worker::reduction(bool i, Depth d, int mn, int delta) const {
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))
    int reductionScale = reductions[d] * reductions[mn];
    return reductionScale - delta * S(731, 731) / rootDelta + !i * reductionScale * S(216, 216) / 512 + S(1089, 1089);
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
Value value_to_tt(Value v, int ply) { return is_win(v) ? v + ply : is_loss(v) ? v - ply : v; }

Value value_from_tt(Value v, int ply, int r50c) {
    if (!is_valid(v)) return VALUE_NONE;
    if (is_win(v)) {
        if (v >= VALUE_MATE_IN_MAX_PLY && VALUE_MATE - v > 100 - r50c) return VALUE_TB_WIN_IN_MAX_PLY - 1;
        if (VALUE_TB - v > 100 - r50c) return VALUE_TB_WIN_IN_MAX_PLY - 1;
        return v - ply;
    }
    if (is_loss(v)) {
        if (v <= VALUE_MATED_IN_MAX_PLY && VALUE_MATE + v > 100 - r50c) return VALUE_TB_LOSS_IN_MAX_PLY + 1;
        if (VALUE_TB + v > 100 - r50c) return VALUE_TB_LOSS_IN_MAX_PLY + 1;
        return v + ply;
    }
    return v;
}

void update_pv(Move* pv, Move move, const Move* childPv) {
    for (*pv++ = move; childPv && *childPv != Move::none();) *pv++ = *childPv++;
    *pv = Move::none();
}

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
        update_quiet_histories<PvNode>(pos, ss, workerThread, bestMove, quietBonus * S(978, 978) / 1024);
        for (Move move : quietsSearched)
            update_quiet_histories<PvNode>(pos, ss, workerThread, move, -quietMalus * S(1115, 1115) / 1024);
    }
    else
    {
        capturedPiece = type_of(pos.piece_on(bestMove.to_sq()));
        captureHistory[movedPiece][bestMove.to_sq()][capturedPiece] << captureBonus * S(1288, 1288) / 1024;
    }

    if (prevSq != SQ_NONE && ((ss - 1)->moveCount == 1 + (ss - 1)->ttHit) && !pos.captured_piece())
        update_continuation_histories<PvNode>(ss - 1, pos.piece_on(prevSq), prevSq, -captureMalus * S(622, 622) / 1024);

    for (Move move : capturesSearched)
    {
        movedPiece    = pos.moved_piece(move);
        capturedPiece = type_of(pos.piece_on(move.to_sq()));
        captureHistory[movedPiece][move.to_sq()][capturedPiece] << -captureMalus * S(1431, 1431) / 1024;
    }
#undef S
}

template<bool PvNode>
void update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus) {
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))
    static constexpr std::array<ConthistBonus, 6> conthist_bonuses = {
      {{1, S(1108, 1108)}, {2, S(652, 652)}, {3, S(273, 273)}, {4, S(572, 572)}, {5, S(126, 126)}, {6, S(449, 449)}}};

    for (const auto [i, weight] : conthist_bonuses)
    {
        if (ss->inCheck && i > 2) break;
        if (((ss - i)->currentMove).is_ok())
            (*(ss - i)->continuationHistory)[pc][to] << (bonus * weight / 1024) + S(80, 80) * (i < 2);
    }
#undef S
}

template<bool PvNode>
void update_quiet_histories(
  const Position& pos, Stack* ss, Search::Worker& workerThread, Move move, int bonus) {
#define S(non_pv, pv) (PvNode ? (pv) : (non_pv))
    Color us = pos.side_to_move();
    workerThread.mainHistory[us][move.from_to()] << bonus;
    if (ss->ply < LOW_PLY_HISTORY_SIZE)
        workerThread.lowPlyHistory[ss->ply][move.from_to()] << (bonus * S(771, 771) / 1024) + S(40, 40);

    update_continuation_histories<PvNode>(ss, pos.moved_piece(move), move.to_sq(),
                                  bonus * (bonus > 0 ? S(979, 979) : S(842, 842)) / 1024);

    int pIndex = pawn_history_index(pos);
    workerThread.pawnHistory[pIndex][pos.moved_piece(move)][move.to_sq()]
      << (bonus * (bonus > 0 ? S(704, 704) : S(439, 439)) / 1024) + S(70, 70);
#undef S
}

}

Move Skill::pick_best(const RootMoves& rootMoves, size_t multiPV) {
    static PRNG rng(now());
    Value  topScore = rootMoves[0].score;
    int    delta    = std::min(topScore - rootMoves[multiPV - 1].score, int(PawnValue));
    int    maxScore = -VALUE_INFINITE;
    double weakness = 120 - 2 * level;

    for (size_t i = 0; i < multiPV; ++i)
    {
        int push = int(weakness * int(topScore - rootMoves[i].score) + delta * (rng.rand<unsigned>() % int(weakness))) / 128;
        if (rootMoves[i].score + push >= maxScore)
        {
            maxScore = rootMoves[i].score + push;
            best     = rootMoves[i].pv[0];
        }
    }
    return best;
}

void SearchManager::check_time(Search::Worker& worker) {
    if (--callsCnt > 0) return;
    callsCnt = worker.limits.nodes ? std::min(512, int(worker.limits.nodes / 1024)) : 512;
    static TimePoint lastInfoTime = now();
    TimePoint elapsed = tm.elapsed([&worker]() { return worker.threads.nodes_searched(); });
    TimePoint tick    = worker.limits.startTime + elapsed;
    if (tick - lastInfoTime >= 1000) { lastInfoTime = tick; dbg_print(); }
    if (ponder) return;
    if (worker.completedDepth >= 1 && ((worker.limits.use_time_management() && (elapsed > tm.maximum() || stopOnPonderhit))
          || (worker.limits.movetime && elapsed >= worker.limits.movetime)
          || (worker.limits.nodes && worker.threads.nodes_searched() >= worker.limits.nodes)))
        worker.threads.stop = worker.threads.abortedSearch = true;
}

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
        return limits.use_time_management() && 2 * std::chrono::duration<double, std::milli>(t_end - t_start).count() > moveOverhead;
    };

    std::list<StateInfo> sts;
    auto& stRoot = sts.emplace_back();
    pos.do_move(rootMove.pv[0], stRoot);
    int ply = 1;

    while (size_t(ply) < rootMove.pv.size()) {
        Move& pvMove = rootMove.pv[ply];
        RootMoves legalMoves;
        for (const auto& m : MoveList<LEGAL>(pos)) legalMoves.emplace_back(m);
        Tablebases::Config config = Tablebases::rank_root_moves(options, pos, legalMoves);
        RootMove&          rm     = *std::find(legalMoves.begin(), legalMoves.end(), pvMove);
        if (legalMoves[0].tbRank != rm.tbRank) break;
        ply++;
        auto& st = sts.emplace_back();
        pos.do_move(pvMove, st);
        if (config.rootInTB && ((rule50 && pos.is_draw(ply)) || pos.is_repetition(ply))) {
            pos.undo_move(pvMove);
            ply--;
            break;
        }
        if (config.rootInTB && time_abort()) break;
    }

    rootMove.pv.resize(ply);

    while (!(rule50 && pos.is_draw(0))) {
        if (time_abort()) break;
        RootMoves legalMoves;
        for (const auto& m : MoveList<LEGAL>(pos)) {
            auto& rm = legalMoves.emplace_back(m);
            StateInfo tmpSI;
            pos.do_move(m, tmpSI);
            for (const auto& mOpp : MoveList<LEGAL>(pos)) rm.tbRank -= pos.capture(mOpp) ? 100 : 1;
            pos.undo_move(m);
        }
        if (legalMoves.size() == 0) break;
        std::stable_sort(legalMoves.begin(), legalMoves.end(), [](const Search::RootMove& a, const Search::RootMove& b) { return a.tbRank > b.tbRank; });
        Tablebases::Config config = Tablebases::rank_root_moves(options, pos, legalMoves, true);
        if (!config.rootInTB || config.cardinality > 0) break;
        ply++;
        Move& pvMove = legalMoves[0].pv[0];
        rootMove.pv.push_back(pvMove);
        auto& st = sts.emplace_back();
        pos.do_move(pvMove, st);
    }

    if (pos.is_draw(0)) v = VALUE_DRAW;
    for (auto it = rootMove.pv.rbegin(); it != rootMove.pv.rend(); ++it) pos.undo_move(*it);
    if (time_abort()) sync_cout << "info string Syzygy based PV extension requires more time, increase Move Overhead as needed." << sync_endl;
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

    for (size_t i = 0; i < multiPV; ++i) {
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
    if (ttHit) {
        if (MoveList<LEGAL>(pos).contains(ttData.move))
            pv.push_back(ttData.move);
    }
    pos.undo_move(pv[0]);
    return pv.size() > 1;
}

}  // namespace Stockfish
