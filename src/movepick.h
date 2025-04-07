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

#ifndef MOVEPICK_H_INCLUDED
#define MOVEPICK_H_INCLUDED

#include "history.h"
#include "movegen.h"
#include "types.h"

namespace Stockfish {

class Position;

class MovePicker {

   public:
    MovePicker(const MovePicker&)            = delete;
    MovePicker& operator=(const MovePicker&) = delete;
    MovePicker(const Position&,
               Move,
               Depth, // Depth helps determine if it's a QSearch call initially
               const ButterflyHistory*,
               const LowPlyHistory*,
               const CapturePieceToHistory*,
               const PieceToHistory**,
               const PawnHistory*,
               int);
    // Constructor for ProbCut (remains unchanged)
    MovePicker(const Position&, Move, int, const CapturePieceToHistory*);

    Move next_move();
    void skip_quiet_moves(); // This might become less relevant for qsearch but keep for main search

   private:
    template<typename Pred>
    Move select(Pred);
    template<GenType>
    void     score();
    ExtMove* begin() { return cur; }
    ExtMove* end() { return endMoves; }

    const Position&              pos;
    const ButterflyHistory*      mainHistory;
    const LowPlyHistory*         lowPlyHistory;
    const CapturePieceToHistory* captureHistory;
    const PieceToHistory**       continuationHistory;
    const PawnHistory*           pawnHistory;
    Move                         ttMove;
    ExtMove *                    cur, *endMoves, *endBadCaptures, *beginBadQuiets, *endBadQuiets;
    int                          stage;
    int                          threshold; // Used by ProbCut constructor
    Depth                        depth;     // Keep depth info if needed elsewhere, e.g., futility
    int                          ply;
    bool                         skipQuiets = false;
    // Flag to indicate if we should only return captures (QSearch behavior)
    bool                         qSearchOnlyCaptures = false; // <-- ADDED
    ExtMove                      moves[MAX_MOVES];
};

}  // namespace Stockfish

#endif  // #ifndef MOVEPICK_H_INCLUDED
