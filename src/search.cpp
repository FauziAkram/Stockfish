#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>

#include "evaluate.h"
#include "misc.h"
#include "movegen.h"
#include "movepick.h"
#include "position.h"
#include "search.h"
#include "thread.h"
#include "timeman.h"
#include "tt.h"
#include "uci.h"
#include "syzygy/tbprobe.h"

namespace Stockfish {

namespace Search {
  LimitsType Limits;
}

namespace Tablebases {
  int Cardinality;
  bool RootInTB;
  bool UseRule50;
  Depth ProbeDepth;
}

namespace TB = Tablebases;
using std::string;
using Eval::evaluate;
using namespace Search;

namespace {

  enum NodeType { NonPV, PV, Root };

  Value futility_margin(Depth d, bool improving) {
    return Value(140 * (d - improving));
  }

  int Reductions[MAX_MOVES];

  Depth reduction(bool i, Depth d, int mn, Value delta, Value rootDelta) {
    int r = Reductions[d] * Reductions[mn];
    return (r + 1372 - int(delta) * 1073 / int(rootDelta)) / 1024 + (!i && r > 936);
  }

  constexpr int futility_move_count(bool improving, Depth depth) {
    return improving ? (3 + depth * depth) : (3 + depth * depth) / 2;
  }

  int stat_bonus(Depth d) {
    return std::min(336 * d - 547, 1561);
  }

  Value value_draw(const Position& pos) {
    return VALUE_DRAW - 1 + Value(pos.this_thread()->nodes & 0x2);
  }

  struct Skill {
    Skill(int skill_level, int uci_elo) {
        if (uci_elo)
        {
            double e = double(uci_elo - 1320) / (3190 - 1320);
            level = std::clamp((((37.2473 * e - 40.8525) * e + 22.2943) * e - 0.311438), 0.0, 19.0);
        }
        else
            level = double(skill_level);
    }
    bool enabled() const { return level < 20.0; }
    bool time_to_pick(Depth depth) const { return depth == 1 + int(level); }
    Move pick_best(size_t multiPV);
    double level;
    Move best = MOVE_NONE;
  };

  template <NodeType nodeType>
  Value search(Position& pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode);

  template <NodeType nodeType>
  Value qsearch(Position& pos, Stack* ss, Value alpha, Value beta, Depth depth = 0);

  Value value_to_tt(Value v, int ply);
  Value value_from_tt(Value v, int ply, int r50c);
  void update_pv(Move* pv, Move move, const Move* childPv);
  void update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus);
  void update_quiet_stats(const Position& pos, Stack* ss, Move move, int bonus);
  void update_all_stats(const Position& pos, Stack* ss, Move bestMove, Value bestValue, Value beta, Square prevSq,
                        Move* quietsSearched, int quietCount, Move* capturesSearched, int captureCount, Depth depth);

  template<bool Root>
  uint64_t perft(Position& pos, Depth depth) {
    StateInfo st;
    uint64_t cnt, nodes = 0;
    const bool leaf = (depth == 2);

    for (const auto& m : MoveList<LEGAL>(pos))
    {
        if (Root && depth <= 1)
            cnt = 1, nodes++;
        else
        {
            pos.do_move(m, st);
            nodes += (cnt = leaf ? MoveList<LEGAL>(pos).size() : perft<false>(pos, depth - 1));
            pos.undo_move(m);
        }
        if (Root)
            sync_cout << UCI::move(m, pos.is_chess960()) << ": " << cnt << sync_endl;
    }
    return nodes;
  }
} // namespace

void Search::init() {
  for (int i = 1; i < MAX_MOVES; ++i)
      Reductions[i] = int((20.57 + std::log(Threads.size()) / 2) * std::log(i));
}

void Search::clear() {
  Threads.main()->wait_for_search_finished();
  Time.availableNodes = 0;
  TT.clear();
  Threads.clear();
  Tablebases::init(Options["SyzygyPath"]);
}

void MainThread::search() {
  if (Limits.perft)
  {
      nodes = perft<true>(rootPos, Limits.perft);
      sync_cout << "\nNodes searched: " << nodes << "\n" << sync_endl;
      return;
  }
  Color us = rootPos.side_to_move();
  Time.init(Limits, us, rootPos.game_ply());
  TT.new_search();

  if (rootMoves.empty())
  {
      rootMoves.emplace_back(MOVE_NONE);
      sync_cout << "info depth 0 score "
                << UCI::value(rootPos.checkers() ? -VALUE_MATE : VALUE_DRAW)
                << sync_endl;
  }
  else
  {
      Threads.start_searching();
      Thread::search();
  }

  while (!Threads.stop && (ponder || Limits.infinite)) {}
  Threads.stop = true;
  Threads.wait_for_search_finished();

  if (Limits.npmsec)
      Time.availableNodes += Limits.inc[us] - Threads.nodes_searched();

  Thread* bestThread = this;
  Skill skill = Skill(Options["Skill Level"], Options["UCI_LimitStrength"] ? int(Options["UCI_Elo"]) : 0);
  if (   int(Options["MultiPV"]) == 1
      && !Limits.depth
      && !skill.enabled()
      && rootMoves[0].pv[0] != MOVE_NONE)
      bestThread = Threads.get_best_thread();

  bestPreviousScore = bestThread->rootMoves[0].score;
  bestPreviousAverageScore = bestThread->rootMoves[0].score; // Simplified

  if (bestThread != this)
      sync_cout << UCI::pv(bestThread->rootPos, bestThread->completedDepth) << sync_endl;

  sync_cout << "bestmove " << UCI::move(bestThread->rootMoves[0].pv[0], rootPos.is_chess960());

  if (bestThread->rootMoves[0].pv.size() > 1 || bestThread->rootMoves[0].extract_ponder_from_tt(rootPos))
      std::cout << " ponder " << UCI::move(bestThread->rootMoves[0].pv[1], rootPos.is_chess960());
  std::cout << sync_endl;
}

void Thread::search() {
  Stack stack[MAX_PLY+10], *ss = stack+7;
  Move  pv[MAX_PLY+1];
  Value alpha, beta, delta;
  Move  lastBestMove = MOVE_NONE;
  Depth lastBestMoveDepth = 0;
  MainThread* mainThread = (this == Threads.main() ? Threads.main() : nullptr);
  double timeReduction = 1, totBestMoveChanges = 0;
  int iterIdx = 0;

  std::memset(ss-7, 0, 10 * sizeof(Stack));
  for (int i = 7; i > 0; --i)
  {
      (ss-i)->continuationHistory = &this->continuationHistory[0][0][NO_PIECE][0];
      (ss-i)->staticEval = VALUE_NONE;
  }
  for (int i = 0; i <= MAX_PLY + 2; ++i)
      (ss+i)->ply = i;
  ss->pv = pv;
  bestValue = -VALUE_INFINITE;

  if (mainThread)
  {
      if (mainThread->bestPreviousScore == VALUE_INFINITE)
          mainThread->iterValue[0] = mainThread->iterValue[1] = mainThread->iterValue[2] = mainThread->iterValue[3] = VALUE_ZERO;
      else
          mainThread->iterValue[0] = mainThread->iterValue[1] = mainThread->iterValue[2] = mainThread->iterValue[3] = mainThread->bestPreviousScore;
  }

  size_t multiPV = size_t(Options["MultiPV"]);
  Skill skill(Options["Skill Level"], Options["UCI_LimitStrength"] ? int(Options["UCI_Elo"]) : 0);
  if (skill.enabled())
      multiPV = std::max(multiPV, (size_t)4);
  multiPV = std::min(multiPV, rootMoves.size());
  int searchAgainCounter = 0;

  while (   ++rootDepth < MAX_PLY
         && !Threads.stop
         && !(Limits.depth && mainThread && rootDepth > Limits.depth))
  {
      if (mainThread)
          totBestMoveChanges /= 2;

      for (RootMove& rm : rootMoves)
          rm.previousScore = rm.score;

      size_t pvFirst = 0;
      pvLast = 0;
      if (!Threads.increaseDepth)
          searchAgainCounter++;

      for (pvIdx = 0; pvIdx < multiPV && !Threads.stop; ++pvIdx)
      {
          if (pvIdx == pvLast)
          {
              pvFirst = pvLast;
              for (pvLast++; pvLast < rootMoves.size(); pvLast++)
                  if (rootMoves[pvLast].tbRank != rootMoves[pvFirst].tbRank)
                      break;
          }
          selDepth = 0;

          Value prev = rootMoves[pvIdx].previousScore;
          delta = Value(10) + int(prev) * prev / 15799;
          alpha = std::max(prev - delta,-VALUE_INFINITE);
          beta  = std::min(prev + delta, VALUE_INFINITE);
          int failedHighCnt = 0;

          while (true)
          {
              Depth adjustedDepth = std::max(1, rootDepth - failedHighCnt - 3 * (searchAgainCounter + 1) / 4);
              bestValue = Stockfish::search<Root>(rootPos, ss, alpha, beta, adjustedDepth, false);
              std::stable_sort(rootMoves.begin() + pvIdx, rootMoves.begin() + pvLast);
              if (Threads.stop)
                  break;
              if (mainThread && multiPV == 1 && (bestValue <= alpha || bestValue >= beta) && Time.elapsed() > 3000)
                  sync_cout << UCI::pv(rootPos, rootDepth) << sync_endl;

              if (bestValue <= alpha)
              {
                  beta = (alpha + beta) / 2;
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
          if (mainThread && (Threads.stop || pvIdx + 1 == multiPV || Time.elapsed() > 3000))
              sync_cout << UCI::pv(rootPos, rootDepth) << sync_endl;
      }
      if (!Threads.stop)
          completedDepth = rootDepth;
      if (rootMoves[0].pv[0] != lastBestMove)
      {
          lastBestMove = rootMoves[0].pv[0];
          lastBestMoveDepth = rootDepth;
      }
      if (Limits.mate && bestValue >= VALUE_MATE_IN_MAX_PLY && VALUE_MATE - bestValue <= 2 * Limits.mate)
          Threads.stop = true;

      if (!mainThread) continue;
      if (skill.enabled() && skill.time_to_pick(rootDepth))
          skill.pick_best(multiPV);

      for (Thread* th : Threads)
      {
          totBestMoveChanges += th->bestMoveChanges;
          th->bestMoveChanges = 0;
      }
      if (Limits.use_time_management() && !Threads.stop && !mainThread->stopOnPonderhit)
      {
          double fallingEval = (69 + 13 * (mainThread->bestPreviousAverageScore - bestValue)
                                    +  6 * (mainThread->iterValue[iterIdx] - bestValue)) / 619.6;
          fallingEval = std::clamp(fallingEval, 0.5, 1.5);
          timeReduction = lastBestMoveDepth + 8 < completedDepth ? 1.57 : 0.65;
          double reduction = (1.4 + mainThread->previousTimeReduction) / (2.08 * timeReduction);
          double bestMoveInstability = 1 + 1.8 * totBestMoveChanges / Threads.size();
          double totalTime = Time.optimum() * fallingEval * reduction * bestMoveInstability;
          if (rootMoves.size() == 1) totalTime = std::min(500.0, totalTime);
          if (Time.elapsed() > totalTime)
          {
              if (mainThread->ponder) mainThread->stopOnPonderhit = true;
              else Threads.stop = true;
          }
          else if (!mainThread->ponder && Time.elapsed() > totalTime * 0.50)
              Threads.increaseDepth = false;
          else
              Threads.increaseDepth = true;
      }
      mainThread->iterValue[iterIdx] = bestValue;
      iterIdx = (iterIdx + 1) & 3;
  }
  if (!mainThread) return;
  mainThread->previousTimeReduction = timeReduction;
  if (skill.enabled())
      std::swap(rootMoves[0], *std::find(rootMoves.begin(), rootMoves.end(), skill.best ? skill.best : skill.pick_best(multiPV)));
}

namespace {
// The rest of the file is the search implementation, adapted from SF16.1
template <NodeType NT>
Value search(Position& pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode) {
  constexpr bool PvNode = NT != NonPV;
  constexpr bool rootNode = NT == Root;
  StateInfo st;
  Move bestMove = MOVE_NONE;
  Thread* thisThread = pos.this_thread();
  
  if (!rootNode)
  {
    if (pos.is_draw(ss->ply) || ss->ply >= MAX_PLY) return value_draw(pos);
    alpha = std::max(mated_in(ss->ply), alpha);
    beta = std::min(mate_in(ss->ply + 1), beta);
    if (alpha >= beta) return alpha;
  }
  if (depth <= 0) return qsearch<NT>(pos, ss, alpha, beta);

  ss->inCheck = pos.checkers();
  (ss + 1)->killers[0] = (ss + 1)->killers[1] = MOVE_NONE;
  (ss + 2)->cutoffCnt = 0;
  ss->doubleExtensions = (ss - 1)->doubleExtensions;

  TTEntry* tte;
  bool ttHit;
  Key posKey = pos.key();
  tte = TT.probe(posKey, ttHit);
  Value ttValue = ttHit ? value_from_tt(tte->value(), ss->ply, pos.rule50_count()) : VALUE_NONE;
  Move ttMove = rootNode ? thisThread->rootMoves[thisThread->pvIdx].pv[0]
                         : (ttHit ? tte->move() : MOVE_NONE);
  ss->ttPv = PvNode || (ttHit && tte->is_pv());
  
  if (!PvNode && ttHit && tte->depth() >= depth && ttValue != VALUE_NONE
      && (tte->bound() & (ttValue >= beta ? BOUND_LOWER : BOUND_UPPER)))
  {
      if (ttMove)
      {
          if (ttValue >= beta)
          {
              if (!pos.capture_stage(ttMove))
                  update_quiet_stats(pos, ss, ttMove, stat_bonus(depth));
              if (to_sq((ss - 1)->currentMove) != SQ_NONE && (ss - 1)->moveCount <= 2 && !pos.captured_piece())
                  update_continuation_histories(ss - 1, pos.piece_on(to_sq((ss-1)->currentMove)), to_sq((ss-1)->currentMove), -stat_bonus(depth + 1));
          }
          else if (!pos.capture_stage(ttMove))
          {
              int penalty = -stat_bonus(depth);
              thisThread->mainHistory[pos.side_to_move()][from_to(ttMove)] << penalty;
              update_continuation_histories(ss, pos.moved_piece(ttMove), to_sq(ttMove), penalty);
          }
      }
      return ttValue;
  }

  if (TB::Cardinality) // Tablebases probe
  {
      int piecesCount = pos.count<ALL_PIECES>();
      if (piecesCount <= TB::Cardinality && (piecesCount < TB::Cardinality || depth >= TB::ProbeDepth) && pos.rule50_count() == 0 && !pos.can_castle(ANY_CASTLING))
      {
          TB::ProbeState err;
          TB::WDLScore wdl = Tablebases::probe_wdl(pos, &err);
          if (err != TB::ProbeState::FAIL)
          {
              thisThread->tbHits.fetch_add(1, std::memory_order_relaxed);
              int drawScore = TB::UseRule50 ? 1 : 0;
              Value value = wdl < -drawScore ? VALUE_MATED_IN_MAX_PLY + ss->ply + 1
                                             : wdl > drawScore ? VALUE_MATE_IN_MAX_PLY - ss->ply - 1
                                                               : VALUE_DRAW + 2 * wdl * drawScore;
              Bound b = wdl < -drawScore ? BOUND_UPPER : wdl > drawScore ? BOUND_LOWER : BOUND_EXACT;
              if (b == BOUND_EXACT || (b == BOUND_LOWER ? value >= beta : value <= alpha))
              {
                  tte->save(posKey, value_to_tt(value, ss->ply), ss->ttPv, b, std::min(MAX_PLY - 1, depth + 6), MOVE_NONE, VALUE_NONE);
                  return value;
              }
              if (PvNode) {
                  if (b == BOUND_LOWER) bestMove = MOVE_NONE, alpha = std::max(alpha, value);
                  // Not updating maxValue here, as TB losses are not searched
              }
          }
      }
  }

  ss->staticEval = ttHit && tte->eval() != VALUE_NONE ? tte->eval() : evaluate(pos);
  if (!ttHit) tte->save(posKey, VALUE_NONE, ss->ttPv, BOUND_NONE, DEPTH_NONE, MOVE_NONE, ss->staticEval);
  Value eval = ss->staticEval;
  if (ttHit && ttValue != VALUE_NONE && (tte->bound() & (ttValue > eval ? BOUND_LOWER : BOUND_UPPER)))
      eval = ttValue;

  bool improving = ss->ply >= 2 && !ss->inCheck && (ss->staticEval > (ss-2)->staticEval || (ss-2)->staticEval == VALUE_NONE);

  if (!PvNode && !ss->inCheck && eval >= beta && eval - futility_margin(depth, improving) >= beta && eval < VALUE_KNOWN_WIN)
      return eval;
  
  if (!PvNode && !ss->inCheck && eval >= beta && pos.non_pawn_material(pos.side_to_move())) {
    Depth R = 3 + depth / 4 + std::min((eval - beta) / 200, 2);
    (ss + 1)->currentMove = MOVE_NULL;
    pos.do_null_move(st);
    Value nullValue = -search<NonPV>(pos, ss + 1, -beta, -beta + 1, depth - R, !cutNode);
    pos.undo_null_move();
    if (nullValue >= beta) {
      if (nullValue >= VALUE_MATE_IN_MAX_PLY) nullValue = beta;
      if (depth < 14) return nullValue;
      Value v = search<NonPV>(pos, ss, beta - 1, beta, depth - R, false);
      if (v >= beta) return nullValue;
    }
  }

  if (PvNode && !ttMove) depth -= 2;
  if (cutNode && depth >= 8 && !ttMove) depth -= 2;

  Value bestValue = -VALUE_INFINITE;
  int moveCount = 0;
  Move capturesSearched[32], quietsSearched[64];
  int captureCount = 0, quietCount = 0;

  const PieceToHistory* contHist[] = { (ss-1)->continuationHistory, (ss-2)->continuationHistory,
                                          nullptr                   , (ss-4)->continuationHistory,
                                          nullptr                   , (ss-6)->continuationHistory };

  Move countermove = (ss-1)->currentMove != MOVE_NONE && to_sq((ss-1)->currentMove) != SQ_NONE ? thisThread->counterMoves[pos.piece_on(to_sq((ss-1)->currentMove))][to_sq((ss-1)->currentMove)] : MOVE_NONE;
  
  MovePicker mp(pos, ttMove, depth, &thisThread->mainHistory, &thisThread->captureHistory, contHist, countermove, ss->killers);

  while (Move move = mp.next_move()) {
      if (!pos.legal(move)) continue;
      if (rootNode && !std::count(thisThread->rootMoves.begin() + thisThread->pvIdx, thisThread->rootMoves.end(), move)) continue;
      
      moveCount++;
      Depth newDepth = depth - 1;
      bool givesCheck = pos.gives_check(move);
      bool is_capture = pos.capture_stage(move);
      ss->currentMove = move;
      ss->continuationHistory = &thisThread->continuationHistory[ss->inCheck][is_capture][pos.moved_piece(move)][to_sq(move)];
      
      // Pruning, LMR, Extensions
      Depth r = 0;
      if (depth >= 2 && moveCount > 1) {
          r = reduction(improving, depth, moveCount, beta - alpha, thisThread->rootDelta);
          if (ss->ttPv) r -= 2;
          if (cutNode) r += 2;
          if (ttMove && pos.capture_stage(ttMove)) r++;
          if ((ss + 1)->cutoffCnt > 3) r++;
          newDepth -= r;
      }
      
      pos.do_move(move, st);
      Value value;
      if (newDepth < 1)
          value = -qsearch<NonPV>(pos, ss + 1, -(alpha + 1), -alpha);
      else if (moveCount > 1)
          value = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode);
      else
          value = -search<NT>(pos, ss + 1, -beta, -alpha, newDepth, false);

      if (r > 0 && value > alpha)
          value = -search<NT>(pos, ss + 1, -beta, -alpha, depth - 1, false);
      
      pos.undo_move(move);
      
      if (Threads.stop) return VALUE_ZERO;
      
      if (value > bestValue) {
          bestValue = value;
          if (value > alpha) {
              bestMove = move;
              alpha = value;
              if (PvNode) update_pv(ss->pv, move, (ss + 1)->pv);
              if (value >= beta) {
                  (ss + 1)->cutoffCnt++;
                  break;
              }
          }
      } else {
          if (is_capture && captureCount < 32) capturesSearched[captureCount++] = move;
          else if (!is_capture && quietCount < 64) quietsSearched[quietCount++] = move;
      }
  }

  if (moveCount == 0) return ss->inCheck ? mated_in(ss->ply) : VALUE_DRAW;
  
  if (bestMove)
      update_all_stats(pos, ss, bestMove, bestValue, beta, to_sq((ss-1)->currentMove), quietsSearched, quietCount, capturesSearched, captureCount, depth);

  tte->save(posKey, value_to_tt(bestValue, ss->ply), ss->ttPv, bestValue >= beta ? BOUND_LOWER : (PvNode && bestMove ? BOUND_EXACT : BOUND_UPPER), depth, bestMove, ss->staticEval);

  return bestValue;
}

// ... qsearch and other helpers from SF16.1 adapted ...
template <NodeType nodeType>
Value qsearch(Position& pos, Stack* ss, Value alpha, Value beta, Depth depth) {
    // Simplified qsearch from SF16
    StateInfo st;
    Value bestValue;
    if (pos.checkers()) {
        ss->inCheck = true;
        bestValue = -VALUE_INFINITE;
    } else {
        ss->inCheck = false;
        bestValue = evaluate(pos);
        if (bestValue >= beta) return bestValue;
        if (bestValue > alpha) alpha = bestValue;
    }

    MovePicker mp(pos, MOVE_NONE, depth, &pos.this_thread()->mainHistory, &pos.this_thread()->captureHistory, nullptr, ss->killers[0]);
    Move move;

    while ((move = mp.next_move(true)) != MOVE_NONE) {
        if (!pos.legal(move)) continue;

        pos.do_move(move, st);
        Value value = -qsearch<nodeType>(pos, ss + 1, -beta, -alpha, depth - 1);
        pos.undo_move(move);

        if (value > bestValue) {
            bestValue = value;
            if (value >= beta) break;
            if (value > alpha) alpha = value;
        }
    }

    if (ss->inCheck && bestValue == -VALUE_INFINITE)
        return mated_in(ss->ply);
    
    return bestValue;
}

Value value_to_tt(Value v, int ply) {
  return v >= VALUE_TB_WIN_IN_MAX_PLY  ? v + ply : v <= VALUE_TB_LOSS_IN_MAX_PLY ? v - ply : v;
}

Value value_from_tt(Value v, int ply, int r50c) {
  if (v == VALUE_NONE) return VALUE_NONE;
  if (v >= VALUE_TB_WIN_IN_MAX_PLY) {
      if (v >= VALUE_MATE_IN_MAX_PLY && VALUE_MATE - v > 99 - r50c) return VALUE_MATE_IN_MAX_PLY - 1;
      return v - ply;
  }
  if (v <= VALUE_TB_LOSS_IN_MAX_PLY) {
      if (v <= VALUE_MATED_IN_MAX_PLY && VALUE_MATE + v > 99 - r50c) return VALUE_MATED_IN_MAX_PLY + 1;
      return v + ply;
  }
  return v;
}

void update_pv(Move* pv, Move move, const Move* childPv) {
  for (*pv++ = move; childPv && *childPv != MOVE_NONE; ) *pv++ = *childPv++;
  *pv = MOVE_NONE;
}

void update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus) {
  for (int i : {1, 2, 4, 6}) {
      if (ss->inCheck && i > 2) break;
      if (is_ok((ss-i)->currentMove)) (*(ss-i)->continuationHistory)[pc][to] << bonus;
  }
}

void update_quiet_stats(const Position& pos, Stack* ss, Move move, int bonus) {
  if (ss->killers[0] != move) {
      ss->killers[1] = ss->killers[0];
      ss->killers[0] = move;
  }
  Thread* thisThread = pos.this_thread();
  thisThread->mainHistory[pos.side_to_move()][from_to(move)] << bonus;
  update_continuation_histories(ss, pos.moved_piece(move), to_sq(move), bonus);
  if (is_ok((ss-1)->currentMove)) {
      thisThread->counterMoves[pos.piece_on(to_sq((ss-1)->currentMove))][to_sq((ss-1)->currentMove)] = move;
  }
}

void update_all_stats(const Position& pos, Stack* ss, Move bestMove, Value bestValue, Value beta, Square prevSq, Move* quietsSearched, int quietCount, Move* capturesSearched, int captureCount, Depth depth) {
    if (!pos.capture_stage(bestMove))
    {
        int bonus = stat_bonus(depth);
        update_quiet_stats(pos, ss, bestMove, bonus);
        for (int i = 0; i < quietCount; ++i) {
            pos.this_thread()->mainHistory[pos.side_to_move()][from_to(quietsSearched[i])] << -bonus;
            update_continuation_histories(ss, pos.moved_piece(quietsSearched[i]), to_sq(quietsSearched[i]), -bonus);
        }
    }
}

Move Skill::pick_best(size_t multiPV) {
  const RootMoves& rootMoves = Threads.main()->rootMoves;
  static PRNG rng(now());
  Value topScore = rootMoves[0].score;
  int delta = std::min(topScore - rootMoves[multiPV - 1].score, PawnValueMg);
  int maxScore = -VALUE_INFINITE;
  double weakness = 120 - 2 * level;
  for (size_t i = 0; i < multiPV; ++i) {
      int push = int((weakness * int(topScore - rootMoves[i].score) + delta * (rng.rand<unsigned>() % int(weakness))) / 128);
      if (rootMoves[i].score + push >= maxScore) {
          maxScore = rootMoves[i].score + push;
          best = rootMoves[i].pv[0];
      }
  }
  return best;
}
} // namespace

void MainThread::check_time() {
  if (--callsCnt > 0) return;
  callsCnt = Limits.nodes ? std::min(1024, int(Limits.nodes / 1024)) : 1024;
  static TimePoint lastInfoTime = now();
  TimePoint elapsed = Time.elapsed();
  TimePoint tick = Limits.startTime + elapsed;
  if (tick - lastInfoTime >= 1000) {
      lastInfoTime = tick;
      dbg_print();
  }
  if (ponder) return;
  if ((Limits.use_time_management() && (elapsed > Time.maximum() - 10 || stopOnPonderhit))
      || (Limits.movetime && elapsed >= Limits.movetime)
      || (Limits.nodes && Threads.nodes_searched() >= (uint64_t)Limits.nodes))
      Threads.stop = true;
}

string UCI::pv(const Position& pos, Depth depth) {
  std::stringstream ss;
  TimePoint elapsed = Time.elapsed() + 1;
  const RootMoves& rootMoves = pos.this_thread()->rootMoves;
  size_t pvIdx = pos.this_thread()->pvIdx;
  size_t multiPV = std::min((size_t)Options["MultiPV"], rootMoves.size());
  uint64_t nodesSearched = Threads.nodes_searched();
  uint64_t tbHits = Threads.tb_hits() + (TB::RootInTB ? rootMoves.size() : 0);

  for (size_t i = 0; i < multiPV; ++i)
  {
      bool updated = rootMoves[i].score != -VALUE_INFINITE;
      if (depth == 1 && !updated && i > 0) continue;
      Depth d = updated ? depth : std::max(1, depth - 1);
      Value v = updated ? rootMoves[i].score : rootMoves[i].previousScore;
      if (v == -VALUE_INFINITE) v = VALUE_ZERO;
      bool tb = TB::RootInTB && abs(v) < VALUE_MATE_IN_MAX_PLY;
      v = tb ? rootMoves[i].tbScore : v;
      if (ss.rdbuf()->in_avail()) ss << "\n";
      ss << "info"
         << " depth "    << d
         << " seldepth " << rootMoves[i].selDepth
         << " multipv "  << i + 1
         << " score "    << UCI::value(v);
      if (Options["UCI_ShowWDL"]) ss << UCI::wdl(v, pos.game_ply());
      ss << " nodes "    << nodesSearched
         << " nps "      << nodesSearched * 1000 / elapsed
         << " hashfull " << TT.hashfull()
         << " tbhits "   << tbHits
         << " time "     << elapsed
         << " pv";
      for (Move m : rootMoves[i].pv) ss << " " << UCI::move(m, pos.is_chess960());
  }
  return ss.str();
}

bool RootMove::extract_ponder_from_tt(Position& pos) {
    StateInfo st;
    bool ttHit;
    assert(pv.size() == 1);
    if (pv[0] == MOVE_NONE) return false;
    pos.do_move(pv[0], st);
    TTEntry* tte = TT.probe(pos.key(), ttHit);
    if (ttHit) {
        Move m = tte->move();
        if (MoveList<LEGAL>(pos).contains(m)) pv.push_back(m);
    }
    pos.undo_move(pv[0]);
    return pv.size() > 1;
}

void Tablebases::rank_root_moves(Position& pos, Search::RootMoves& rootMoves) {
    RootInTB = false;
    UseRule50 = bool(Options["Syzygy50MoveRule"]);
    ProbeDepth = int(Options["SyzygyProbeDepth"]);
    Cardinality = int(Options["SyzygyProbeLimit"]);
    bool dtz_available = true;
    if (Cardinality > MaxCardinality) {
        Cardinality = MaxCardinality;
        ProbeDepth = 0;
    }
    if (Cardinality >= popcount(pos.pieces()) && !pos.can_castle(ANY_CASTLING)) {
        RootInTB = root_probe(pos, rootMoves);
        if (!RootInTB) {
            dtz_available = false;
            RootInTB = root_probe_wdl(pos, rootMoves);
        }
    }
    if (RootInTB) {
        std::stable_sort(rootMoves.begin(), rootMoves.end(), [](const RootMove &a, const RootMove &b) { return a.tbRank > b.tbRank; } );
        if (dtz_available || rootMoves[0].tbScore <= VALUE_DRAW) Cardinality = 0;
    } else {
        for (auto& m : rootMoves) m.tbRank = 0;
    }
}

} // namespace Stockfish
