#ifndef MOVEPICK_H_INCLUDED
#define MOVEPICK_H_INCLUDED

#include <array>
#include <limits>
#include <type_traits>

#include "movegen.h"
#include "position.h"
#include "types.h"

namespace Stockfish {

// The StatsEntry/Stats templates from SF16.1 are what we'll use for history tables.
template<typename T, int D>
class StatsEntry {
  T entry;
public:
  void operator=(const T& v) { entry = v; }
  T* operator&() { return &entry; }
  T* operator->() { return &entry; }
  operator const T&() const { return entry; }
  void operator<<(int bonus) {
    assert(abs(bonus) <= D);
    static_assert(D <= std::numeric_limits<T>::max(), "D overflows T");
    entry += bonus - entry * abs(bonus) / D;
    assert(abs(entry) <= D);
  }
};

template <typename T, int D, int Size, int... Sizes>
struct Stats : public std::array<Stats<T, D, Sizes...>, Size> {
  using stats = Stats<T, D, Size, Sizes...>;
  void fill(const T& v) {
    assert(std::is_standard_layout<stats>::value);
    using entry = StatsEntry<T, D>;
    entry* p = reinterpret_cast<entry*>(this);
    std::fill(p, p + sizeof(*this) / sizeof(entry), v);
  }
};

template <typename T, int D, int Size>
struct Stats<T, D, Size> : public std::array<StatsEntry<T, D>, Size> {};

enum StatsParams { NOT_USED = 0 };
enum StatsType { NoCaptures, Captures };

using ButterflyHistory = Stats<int16_t, 7183, COLOR_NB, int(SQUARE_NB) * int(SQUARE_NB)>;
using CounterMoveHistory = Stats<Move, NOT_USED, PIECE_NB, SQUARE_NB>;
using CapturePieceToHistory = Stats<int16_t, 10692, PIECE_NB, SQUARE_NB, PIECE_TYPE_NB>;
using PieceToHistory = Stats<int16_t, 29952, PIECE_NB, SQUARE_NB>;
using ContinuationHistory = Stats<PieceToHistory, NOT_USED, PIECE_NB, SQUARE_NB>;

class MovePicker {
  enum PickType { Next, Best };
public:
  MovePicker(const MovePicker&) = delete;
  MovePicker& operator=(const MovePicker&) = delete;
  MovePicker(const Position&, Move, Depth, const ButterflyHistory*,
                                           const CapturePieceToHistory*,
                                           const PieceToHistory**,
                                           Move,
                                           const Move*);
  MovePicker(const Position&, Move, Depth, const ButterflyHistory*,
                                           const CapturePieceToHistory*,
                                           const PieceToHistory**,
                                           Square);
  MovePicker(const Position&, Move, Value, const CapturePieceToHistory*);
  Move next_move(bool skipQuiets = false);

private:
  template<PickType T, typename Pred> Move select(Pred);
  template<GenType> void score();
  ExtMove* begin() { return cur; }
  ExtMove* end() { return endMoves; }

  const Position& pos;
  const ButterflyHistory* mainHistory;
  const CapturePieceToHistory* captureHistory;
  const PieceToHistory** continuationHistory;
  Move ttMove;
  ExtMove refutations[3], *cur, *endMoves, *endBadCaptures;
  int stage;
  Square recaptureSquare;
  Value threshold;
  Depth depth;
  ExtMove moves[MAX_MOVES];
};

} // namespace Stockfish

#endif // #ifndef MOVEPICK_H_INCLUDED
