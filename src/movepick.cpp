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

#include "movepick.h"

#include <cassert>
#include <limits> // Required for numeric_limits

#include "bitboard.h"
#include "misc.h"
#include "position.h"

namespace Stockfish {

namespace {

enum Stages {
    // generate main search moves
    MAIN_TT,
    CAPTURE_INIT,
    GOOD_CAPTURE,
    QUIET_INIT,
    GOOD_QUIET,
    BAD_CAPTURE,
    BAD_QUIET,

    // generate evasion moves
    EVASION_TT,
    EVASION_INIT,
    EVASION,

    // generate probcut moves
    PROBCUT_TT,
    PROBCUT_INIT,
    PROBCUT

    // Removed QSearch stages:
    // QSEARCH_TT,
    // QCAPTURE_INIT,
    // QCAPTURE
};

// Sort moves in descending order up to and including a given limit.
// The order of moves smaller than the limit is left unspecified.
void partial_insertion_sort(ExtMove* begin, ExtMove* end, int limit) {
    // (Implementation remains the same)
    for (ExtMove *sortedEnd = begin, *p = begin + 1; p < end; ++p)
        if (p->value >= limit)
        {
            ExtMove tmp = *p, *q;
            *p          = *++sortedEnd;
            for (q = sortedEnd; q != begin && *(q - 1) < tmp; --q)
                *q = *(q - 1);
            *q = tmp;
        }
}

} // namespace


// MovePicker constructor for the main search and (now also) for the quiescence search
MovePicker::MovePicker(const Position&              p,
                       Move                         ttm,
                       Depth                        d,
                       const ButterflyHistory*      mh,
                       const LowPlyHistory*         lph,
                       const CapturePieceToHistory* cph,
                       const PieceToHistory**       ch,
                       const PawnHistory*           ph,
                       int                          pl) :
    pos(p),
    mainHistory(mh),
    lowPlyHistory(lph),
    captureHistory(cph),
    continuationHistory(ch),
    pawnHistory(ph),
    ttMove(ttm),
    depth(d), // Keep depth, might be used for thresholds etc.
    ply(pl) {

    assert(!pos.checkers() || d > 0); // QSearch should not be called when in check

    if (pos.checkers()) {
        // Evasions path remains the same
        stage = EVASION_TT + !(ttm && pos.pseudo_legal(ttm));
    } else {
        // Always start with main search stages now
        stage = MAIN_TT + !(ttm && pos.pseudo_legal(ttm));
        // Set the flag if this constructor is called in a QSearch context
        if (depth <= 0) {
            qSearchOnlyCaptures = true;
        }
    }
}

// MovePicker constructor for ProbCut (remains unchanged)
MovePicker::MovePicker(const Position& p, Move ttm, int th, const CapturePieceToHistory* cph) :
    pos(p),
    captureHistory(cph),
    ttMove(ttm),
    threshold(th) {
    assert(!pos.checkers());
    stage = PROBCUT_TT
          + !(ttm && pos.capture_stage(ttm) && pos.pseudo_legal(ttm) && pos.see_ge(ttm, threshold));
}


// Scoring function remains unchanged, it's generic based on GenType
template<GenType Type>
void MovePicker::score() {
    // (Implementation remains the same)
    static_assert(Type == CAPTURES || Type == QUIETS || Type == EVASIONS, "Wrong type");

    [[maybe_unused]] Bitboard threatenedByPawn, threatenedByMinor, threatenedByRook,
      threatenedPieces;
    if constexpr (Type == QUIETS)
    {
        Color us = pos.side_to_move();

        threatenedByPawn = pos.attacks_by<PAWN>(~us);
        threatenedByMinor =
          pos.attacks_by<KNIGHT>(~us) | pos.attacks_by<BISHOP>(~us) | threatenedByPawn;
        threatenedByRook = pos.attacks_by<ROOK>(~us) | threatenedByMinor;

        // Pieces threatened by pieces of lesser material value
        threatenedPieces = (pos.pieces(us, QUEEN) & threatenedByRook)
                         | (pos.pieces(us, ROOK) & threatenedByMinor)
                         | (pos.pieces(us, KNIGHT, BISHOP) & threatenedByPawn);
    }

    for (auto& m : *this)
        if constexpr (Type == CAPTURES)
            m.value =
              7 * int(PieceValue[pos.piece_on(m.to_sq())])
              + (*captureHistory)[pos.moved_piece(m)][m.to_sq()][type_of(pos.piece_on(m.to_sq()))];

        else if constexpr (Type == QUIETS)
        {
            Piece     pc   = pos.moved_piece(m);
            PieceType pt   = type_of(pc);
            Square    from = m.from_sq();
            Square    to   = m.to_sq();

            // histories
            m.value = 2 * (*mainHistory)[pos.side_to_move()][m.from_to()];
            m.value += 2 * (*pawnHistory)[pawn_structure_index(pos)][pc][to];
            m.value += (*continuationHistory[0])[pc][to];
            m.value += (*continuationHistory[1])[pc][to];
            m.value += (*continuationHistory[2])[pc][to];
            m.value += (*continuationHistory[3])[pc][to];
            m.value += (*continuationHistory[4])[pc][to] / 3;
            m.value += (*continuationHistory[5])[pc][to];

            // bonus for checks
            m.value += bool(pos.check_squares(pt) & to) * 16384;

            // bonus for escaping from capture
            m.value += threatenedPieces & from ? (pt == QUEEN && !(to & threatenedByRook)   ? 51700
                                                  : pt == ROOK && !(to & threatenedByMinor) ? 25600
                                                  : !(to & threatenedByPawn)                ? 14450
                                                                                            : 0)
                                               : 0;

            // malus for putting piece en prise
            m.value -= (pt == QUEEN && bool(to & threatenedByRook)   ? 49000
                        : pt == ROOK && bool(to & threatenedByMinor) ? 24335
                                                                     : 0);

            if (ply < LOW_PLY_HISTORY_SIZE)
                m.value += 8 * (*lowPlyHistory)[ply][m.from_to()] / (1 + 2 * ply);
        }

        else  // Type == EVASIONS
        {
            if (pos.capture_stage(m))
                m.value = PieceValue[pos.piece_on(m.to_sq())] + (1 << 28);
            else
                m.value = (*mainHistory)[pos.side_to_move()][m.from_to()]
                        + (*continuationHistory[0])[pos.moved_piece(m)][m.to_sq()]
                        + (*pawnHistory)[pawn_structure_index(pos)][pos.moved_piece(m)][m.to_sq()];
        }
}

// Select function remains unchanged
template<typename Pred>
Move MovePicker::select(Pred filter) {
    // (Implementation remains the same)
    for (; cur < endMoves; ++cur)
        if (*cur != ttMove && filter())
            return *cur++;
    return Move::none();
}


Move MovePicker::next_move() {

    auto quiet_threshold = [](Depth d) { return -3560 * d; };

top: // Label for goto jumps
    switch (stage)
    {

    case MAIN_TT : // Used for both main search and QSearch now
    case EVASION_TT :
    // Removed QSEARCH_TT case
    case PROBCUT_TT :
        ++stage;
        return ttMove; // Return TT move first if valid

    case CAPTURE_INIT :
    case PROBCUT_INIT :
    // Removed QCAPTURE_INIT case
        cur = endBadCaptures = moves;
        endMoves             = generate<CAPTURES>(pos, cur);
        score<CAPTURES>();
        // Sort all captures initially; good/bad split happens later
        partial_insertion_sort(cur, endMoves, std::numeric_limits<int>::min());
        // Set next stage based on context (ProbCut or regular/QSearch)
        stage = (stage == PROBCUT_INIT) ? PROBCUT : GOOD_CAPTURE;
        goto top; // Re-evaluate switch with the new stage

    case GOOD_CAPTURE :
        if (select([&]() {
                // Use SEE to filter good captures
                return pos.see_ge(*cur, -cur->value / 18) ? true
                                                          : (*endBadCaptures++ = *cur, false);
            }))
            return *(cur - 1); // Return the selected good capture

        // End of good captures. Decide whether to proceed to quiets or bad captures.
        if (qSearchOnlyCaptures) {
             // If QSearch context, skip quiets and go directly to bad captures
            stage = BAD_CAPTURE;
            cur = moves; // Reset pointers to the beginning of the buffer for bad captures
            endMoves = endBadCaptures; // Set end pointer for bad captures
        } else {
            // If main search context, proceed to initialize quiet moves
            stage = QUIET_INIT;
        }
        goto top; // Re-evaluate switch with the new stage

    case QUIET_INIT :
        // This stage is skipped if qSearchOnlyCaptures is true
        assert(!qSearchOnlyCaptures);
        if (!skipQuiets)
        {
            cur = endBadCaptures; // Start generating quiets after the bad captures area
            endMoves = beginBadQuiets = endBadQuiets = generate<QUIETS>(pos, cur);
            score<QUIETS>();
            partial_insertion_sort(cur, endMoves, quiet_threshold(depth));
        }
        stage = GOOD_QUIET; // Proceed to good quiets
        goto top; // Re-evaluate switch

    case GOOD_QUIET :
        // This stage is skipped if qSearchOnlyCaptures is true
        assert(!qSearchOnlyCaptures);
        if (!skipQuiets && select([]() { return true; }))
        {
            // Check if the selected quiet move is actually bad based on threshold
             if ((cur - 1)->value > -7998 || (cur - 1)->value <= quiet_threshold(depth))
                return *(cur - 1); // Return good quiet

            // Remaining quiets are considered bad
            beginBadQuiets = cur - 1;
        }
        // Prepare the pointers to loop over the bad captures next
        cur = moves;
        endMoves = endBadCaptures;
        stage = BAD_CAPTURE; // Proceed to bad captures
        goto top; // Re-evaluate switch

    case BAD_CAPTURE :
        if (select([]() { return true; })) // Select next bad capture
            return *(cur - 1);

        // End of bad captures. Decide whether to proceed to bad quiets or finish.
        if (qSearchOnlyCaptures) {
            // If QSearch context, we are done after bad captures
            return Move::none();
        } else {
             // If main search context, proceed to bad quiets
            cur = beginBadQuiets; // Prepare pointers for bad quiets
            endMoves = endBadQuiets;
            stage = BAD_QUIET;
            goto top; // Re-evaluate switch
        }

    case BAD_QUIET :
        // This stage is skipped if qSearchOnlyCaptures is true
        assert(!qSearchOnlyCaptures);
        if (!skipQuiets)
            return select([]() { return true; }); // Select next bad quiet

        // End of all moves for main search
        return Move::none();

    // --- Evasion Stages (Unchanged) ---
    case EVASION_INIT :
        cur      = moves;
        endMoves = generate<EVASIONS>(pos, cur);
        score<EVASIONS>();
        partial_insertion_sort(cur, endMoves, std::numeric_limits<int>::min());
        stage = EVASION;
        goto top;

    case EVASION :
        return select([]() { return true; });

    // --- ProbCut Stage (Unchanged) ---
    case PROBCUT :
        return select([&]() { return pos.see_ge(*cur, threshold); });

    // Removed QCAPTURE case
    }

    assert(false); // Should not be reachable
    return Move::none();
}

// skip_quiet_moves function remains unchanged
void MovePicker::skip_quiet_moves() { skipQuiets = true; }

} // namespace Stockfish
