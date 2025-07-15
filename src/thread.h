#ifndef THREAD_H_INCLUDED
#define THREAD_H_INCLUDED

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

#include "material.h"
#include "movepick.h"
#include "pawns.h"
#include "position.h"
#include "search.h"
#include "thread_win32_osx.h"

namespace Stockfish {

class Thread {

  std::mutex mutex;
  std::condition_variable cv;
  size_t idx;
  bool exit = false, searching = true;
  NativeThread stdThread;

public:
  explicit Thread(size_t);
  virtual ~Thread();
  virtual void search();
  void clear();
  void idle_loop();
  void start_searching();
  void wait_for_search_finished();
  size_t id() const { return idx; }

  // HCE-related history tables
  Pawns::Table pawnsTable;
  Material::Table materialTable;
  CounterMoveHistory counterMoves;
  ButterflyHistory mainHistory;
  CapturePieceToHistory captureHistory;
  ContinuationHistory continuationHistory[2][2];

  size_t pvIdx, pvLast;
  std::atomic<uint64_t> nodes, tbHits, bestMoveChanges;
  int selDepth, nmpMinPly;
  Value bestValue; // Used by HCE lazy evaluation

  Position rootPos;
  StateInfo rootState;
  Search::RootMoves rootMoves;
  Depth rootDepth, completedDepth;
  Value rootDelta;
};

struct MainThread : public Thread {
  using Thread::Thread;
  void search() override;
  void check_time();

  double previousTimeReduction;
  Value bestPreviousScore;
  Value bestPreviousAverageScore;
  Value iterValue[4];
  int callsCnt;
  bool stopOnPonderhit;
  std::atomic_bool ponder;
};

struct ThreadPool {
  void start_thinking(Position&, StateListPtr&, const Search::LimitsType&, bool = false);
  void clear();
  void set(size_t);

  MainThread* main()        const { return static_cast<MainThread*>(threads.front()); }
  uint64_t nodes_searched() const { return accumulate(&Thread::nodes); }
  uint64_t tb_hits()        const { return accumulate(&Thread::tbHits); }
  Thread* get_best_thread() const;
  void start_searching();
  void wait_for_search_finished() const;

  std::atomic_bool stop, increaseDepth;

  auto cbegin() const noexcept { return threads.cbegin(); }
  auto begin() noexcept { return threads.begin(); }
  auto end() noexcept { return threads.end(); }
  auto cend() const noexcept { return threads.cend(); }
  auto size() const noexcept { return threads.size(); }
  auto empty() const noexcept { return threads.empty(); }

private:
  StateListPtr setupStates;
  std::vector<Thread*> threads;

  uint64_t accumulate(std::atomic<uint64_t> Thread::* member) const {
    uint64_t sum = 0;
    for (Thread* th : threads)
        sum += (th->*member).load(std::memory_order_relaxed);
    return sum;
  }
};

extern ThreadPool Threads;

} // namespace Stockfish

#endif // #ifndef THREAD_H_INCLUDED
