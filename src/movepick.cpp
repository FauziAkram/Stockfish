#include <cassert>

#include "bitboard.h"
#include "movepick.h"
#include "position.h" // Add this include

namespace Stockfish {

namespace {

  enum Stages {
    MAIN_TT, CAPTURE_INIT, GOOD_CAPTURE, REFUTATION, QUIET_INIT, QUIET, BAD_CAPTURE,
    EVASION_TT, EVASION_INIT, EVASION,
    PROBCUT_TT, PROBCUT_INIT, PROBCUT,
    QSEARCH_TT, QCAPTURE_INIT, QCAPTURE, QCHECK_INIT, QCHECK
  };

  void partial_insertion_sort(ExtMove* begin, ExtMove* end, int limit) {
    for (ExtMove *sortedEnd = begin, *p = begin + 1; p < end; ++p)
        if (p->value >= limit)
        {
            ExtMove tmp = *p, *q;
            *p = *++sortedEnd;
            for (q = sortedEnd; q != begin && *(q - 1) < tmp; --q)
                *q = *(q - 1);
            *q = tmp;
        }
  }

} // namespace

MovePicker::MovePicker(const Position& p, Move ttm, Depth d, const ButterflyHistory* mh,
                                                             const CapturePieceToHistory* cph,
                                                             const PieceToHistory** ch,
                                                             Move cm,
                                                             const Move* killers)
           : pos(p), mainHistory(mh), captureHistory(cph), continuationHistory(ch),
             ttMove(ttm), refutations{{killers[0], 0}, {killers[1], 0}, {cm, 0}}, depth(d)
{
  assert(d > 0);
  stage = (pos.checkers() ? EVASION_TT : MAIN_TT) +
          !(ttm && pos.pseudo_legal(ttm));
}

MovePicker::MovePicker(const Position& p, Move ttm, Depth d, const ButterflyHistory* mh,
                                                             const CapturePieceToHistory* cph,
                                                             const PieceToHistory** ch,
                                                             Square rs)
           : pos(p), mainHistory(mh), captureHistory(cph), continuationHistory(ch), ttMove(ttm), recaptureSquare(rs), depth(d)
{
  assert(d <= 0);
  stage = (pos.checkers() ? EVASION_TT : QSEARCH_TT) +
          !(   ttm
            && pos.pseudo_legal(ttm));
}

MovePicker::MovePicker(const Position& p, Move ttm, Value th, const CapturePieceToHistory* cph)
           : pos(p), captureHistory(cph), ttMove(ttm), threshold(th)
{
  assert(!pos.checkers());
  stage = PROBCUT_TT + !(ttm && pos.capture_stage(ttm)
                             && pos.pseudo_legal(ttm)
                             && pos.see_ge(ttm, threshold));
}

template<GenType Type>
void MovePicker::score() {
  static_assert(Type == CAPTURES || Type == QUIETS || Type == EVASIONS, "Wrong type");
  [[maybe_unused]] Bitboard threatenedByPawn, threatenedByMinor, threatenedByRook, threatenedPieces;
  if constexpr (Type == QUIETS)
  {
      Color us = pos.side_to_move();
      threatenedByPawn  = pos.attacks_by<PAWN>(~us);
      threatenedByMinor = pos.attacks_by<KNIGHT>(~us) | pos.attacks_by<BISHOP>(~us) | threatenedByPawn;
      threatenedByRook  = pos.attacks_by<ROOK>(~us) | threatenedByMinor;
      threatenedPieces = (pos.pieces(us, QUEEN) & threatenedByRook)
                       | (pos.pieces(us, ROOK)  & threatenedByMinor)
                       | (pos.pieces(us, KNIGHT, BISHOP) & threatenedByPawn);
  }
  for (auto& m : *this)
      if constexpr (Type == CAPTURES)
          m.value =  (7 * int(PieceValue[MG][pos.piece_on(to_sq(m))])
                   +     (*captureHistory)[pos.moved_piece(m)][to_sq(m)][type_of(pos.piece_on(to_sq(m)))]) / 16;
      else if constexpr (Type == QUIETS)
          m.value =  2 * (*mainHistory)[pos.side_to_move()][from_to(m)]
                   + 2 * (*continuationHistory[0])[pos.moved_piece(m)][to_sq(m)]
                   +     (*continuationHistory[1])[pos.moved_piece(m)][to_sq(m)]
                   +     (*continuationHistory[3])[pos.moved_piece(m)][to_sq(m)]
                   +     (*continuationHistory[5])[pos.moved_piece(m)][to_sq(m)]
                   +     (threatenedPieces & from_sq(m) ?
                           (type_of(pos.moved_piece(m)) == QUEEN && !(to_sq(m) & threatenedByRook)  ? 50000
                          : type_of(pos.moved_piece(m)) == ROOK  && !(to_sq(m) & threatenedByMinor) ? 25000
                          :                                         !(to_sq(m) & threatenedByPawn)  ? 15000
                          :                                                                           0)
                          :                                                                           0)
                   +     bool(pos.check_squares(type_of(pos.moved_piece(m))) & to_sq(m)) * 16384;
      else // Type == EVASIONS
      {
          if (pos.capture_stage(m))
              m.value =  PieceValue[MG][pos.piece_on(to_sq(m))]
                       - Value(type_of(pos.moved_piece(m)))
                       + (1 << 28);
          else
              m.value =  (*mainHistory)[pos.side_to_move()][from_to(m)]
                       + (*continuationHistory[0])[pos.moved_piece(m)][to_sq(m)];
      }
}

template<MovePicker::PickType T, typename Pred>
Move MovePicker::select(Pred filter) {
  while (cur < endMoves)
  {
      if constexpr (T == Best)
          std::swap(*cur, *std::max_element(cur, endMoves));
      if (*cur != ttMove && filter())
          return *cur++;
      cur++;
  }
  return MOVE_NONE;
}

Move MovePicker::next_move(bool skipQuiets) {
top:
  switch (stage) {
  case MAIN_TT:
  case EVASION_TT:
  case QSEARCH_TT:
  case PROBCUT_TT:
      ++stage;
      return ttMove;
  case CAPTURE_INIT:
  case PROBCUT_INIT:
  case QCAPTURE_INIT:
      cur = endBadCaptures = moves;
      endMoves = generate<CAPTURES>(pos, cur);
      score<CAPTURES>();
      partial_insertion_sort(cur, endMoves, std::numeric_limits<int>::min());
      ++stage;
      goto top;
  case GOOD_CAPTURE:
      if (select<Next>([&](){
                       return pos.see_ge(*cur, -cur->value) ?
                              true : (*endBadCaptures++ = *cur, false); }))
          return *(cur - 1);
      cur = std::begin(refutations);
      endMoves = std::end(refutations);
      if (   refutations[0].move == refutations[2].move
          || refutations[1].move == refutations[2].move)
          --endMoves;
      ++stage;
      [[fallthrough]];
  case REFUTATION:
      if (select<Next>([&](){ return    *cur != MOVE_NONE
                                    && !pos.capture_stage(*cur)
                                    &&  pos.pseudo_legal(*cur); }))
          return *(cur - 1);
      ++stage;
      [[fallthrough]];
  case QUIET_INIT:
      if (!skipQuiets)
      {
          cur = endBadCaptures;
          endMoves = generate<QUIETS>(pos, cur);
          score<QUIETS>();
          partial_insertion_sort(cur, endMoves, -3000 * depth);
      }
      ++stage;
      [[fallthrough]];
  case QUIET:
      if (   !skipQuiets
          && select<Next>([&](){return   *cur != refutations[0].move
                                      && *cur != refutations[1].move
                                      && *cur != refutations[2].move;}))
          return *(cur - 1);
      cur = moves;
      endMoves = endBadCaptures;
      ++stage;
      [[fallthrough]];
  case BAD_CAPTURE:
      return select<Next>([](){ return true; });
  case EVASION_INIT:
      cur = moves;
      endMoves = generate<EVASIONS>(pos, cur);
      score<EVASIONS>();
      ++stage;
      [[fallthrough]];
  case EVASION:
      return select<Best>([](){ return true; });
  case PROBCUT:
      return select<Next>([&](){ return pos.see_ge(*cur, threshold); });
  case QCAPTURE:
      if (select<Next>([&](){ return   depth > DEPTH_QS_RECAPTURES
                                    || to_sq(*cur) == recaptureSquare; }))
          return *(cur - 1);
      if (depth != DEPTH_QS_CHECKS)
          return MOVE_NONE;
      ++stage;
      [[fallthrough]];
  case QCHECK_INIT:
      cur = moves;
      endMoves = generate<QUIET_CHECKS>(pos, cur);
      ++stage;
      [[fallthrough]];
  case QCHECK:
      return select<Next>([](){ return true; });
  }
  assert(false);
  return MOVE_NONE;
}

} // namespace Stockfish
