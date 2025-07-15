#include <cassert>

#include <algorithm> // For std::count
#include "movegen.h"
#include "search.h"
#include "thread.h"
#include "uci.h"
#include "syzygy/tbprobe.h"
#include "tt.h"

namespace Stockfish {

ThreadPool Threads;

Thread::Thread(size_t n) : idx(n), stdThread(&Thread::idle_loop, this) {
  wait_for_search_finished();
}

Thread::~Thread() {
  assert(!searching);
  exit = true;
  start_searching();
  stdThread.join();
}

void Thread::clear() {
  counterMoves.fill(MOVE_NONE);
  mainHistory.fill(0);
  captureHistory.fill(0);

  for (bool inCheck : { false, true })
      for (StatsType c : { NoCaptures, Captures })
          for (auto& to : continuationHistory[inCheck][c])
              for (auto& h : to)
                  h->fill(-71);
}

void Thread::start_searching() {
  mutex.lock();
  searching = true;
  mutex.unlock();
  cv.notify_one();
}

void Thread::wait_for_search_finished() {
  std::unique_lock<std::mutex> lk(mutex);
  cv.wait(lk, [&]{ return !searching; });
}

void Thread::idle_loop() {
  if (Options["Threads"] > 8)
      WinProcGroup::bindThisThread(idx);
  while (true)
  {
      std::unique_lock<std::mutex> lk(mutex);
      searching = false;
      cv.notify_one();
      cv.wait(lk, [&]{ return searching; });
      if (exit)
          return;
      lk.unlock();
      search();
  }
}

void ThreadPool::set(size_t requested) {
  if (threads.size() > 0)
  {
      main()->wait_for_search_finished();
      while (threads.size() > 0)
          delete threads.back(), threads.pop_back();
  }
  if (requested > 0)
  {
      threads.push_back(new MainThread(0));
      while (threads.size() < requested)
          threads.push_back(new Thread(threads.size()));
      clear();
      TT.resize(size_t(Options["Hash"]));
      Search::init();
  }
}

void ThreadPool::clear() {
  for (Thread* th : threads)
      th->clear();
  main()->callsCnt = 0;
  main()->bestPreviousScore = VALUE_INFINITE;
  main()->bestPreviousAverageScore = VALUE_INFINITE;
  main()->previousTimeReduction = 1.0;
}

void ThreadPool::start_thinking(Position& pos, StateListPtr& states,
                                const Search::LimitsType& limits, bool ponderMode) {
  main()->wait_for_search_finished();
  main()->stopOnPonderhit = stop = false;
  increaseDepth = true;
  main()->ponder = ponderMode;
  Search::Limits = limits;
  Search::RootMoves rootMoves;
  for (const auto& m : MoveList<LEGAL>(pos))
      if (   limits.searchmoves.empty()
          || std::count(limits.searchmoves.begin(), limits.searchmoves.end(), m))
          rootMoves.emplace_back(m);
  if (!rootMoves.empty())
      Tablebases::rank_root_moves(pos, rootMoves);
  assert(states.get() || setupStates.get());
  if (states.get())
      setupStates = std::move(states);
  for (Thread* th : threads)
  {
      th->nodes = th->tbHits = th->nmpMinPly = th->bestMoveChanges = 0;
      th->rootDepth = th->completedDepth = 0;
      th->rootMoves = rootMoves;
      th->rootPos.set(pos.fen(), pos.is_chess960(), &th->rootState);
      th->rootState = setupStates->back();
      // Add a back-link from position to thread
      th->rootPos.thisThread = th;
  }
  main()->start_searching();
}

Thread* ThreadPool::get_best_thread() const {
    Thread* bestThread = threads.front();
    std::map<Move, int64_t> votes;
    Value minScore = VALUE_NONE;

    for (Thread* th: threads)
        minScore = std::min(minScore, th->rootMoves[0].score);

    auto thread_value = [minScore](Thread* th) {
            return (th->rootMoves[0].score - minScore + 14) * int(th->completedDepth);
        };

    for (Thread* th : threads)
        votes[th->rootMoves[0].pv[0]] += thread_value(th);

    for (Thread* th : threads)
        if (abs(bestThread->rootMoves[0].score) >= VALUE_TB_WIN_IN_MAX_PLY)
        {
            if (th->rootMoves[0].score > bestThread->rootMoves[0].score)
                bestThread = th;
        }
        else if (   th->rootMoves[0].score >= VALUE_TB_WIN_IN_MAX_PLY
                 || (   th->rootMoves[0].score > VALUE_TB_LOSS_IN_MAX_PLY
                     && (   votes[th->rootMoves[0].pv[0]] > votes[bestThread->rootMoves[0].pv[0]]
                         || (   votes[th->rootMoves[0].pv[0]] == votes[bestThread->rootMoves[0].pv[0]]
                             &&   thread_value(th) * int(th->rootMoves[0].pv.size() > 2)
                                > thread_value(bestThread) * int(bestThread->rootMoves[0].pv.size() > 2)))))
            bestThread = th;
    return bestThread;
}

void ThreadPool::start_searching() {
    for (Thread* th : threads)
        if (th != threads.front())
            th->start_searching();
}

void ThreadPool::wait_for_search_finished() const {
    for (Thread* th : threads)
        if (th != threads.front())
            th->wait_for_search_finished();
}

} // namespace Stockfish
