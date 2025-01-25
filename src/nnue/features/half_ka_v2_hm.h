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

//Definition of input features HalfKP of NNUE evaluation function

#ifndef NNUE_FEATURES_HALF_KA_V2_HM_H_INCLUDED
#define NNUE_FEATURES_HALF_KA_V2_HM_H_INCLUDED

#include <cstdint>

#include "../../misc.h"
#include "../../types.h"
#include "../nnue_common.h"

namespace Stockfish {
struct StateInfo;
class Position;
}

namespace Stockfish::Eval::NNUE::Features {

// Feature HalfKAv2_hm: Combination of the position of own king and the
// position of pieces. Position mirrored such that king is always on e..h files.
class HalfKAv2_hm {

    // Unique number for each piece type on each square
    enum {
        PS_NONE     = 0,
        PS_W_PAWN   = 0,
        PS_B_PAWN   = 1 * SQUARE_NB,
        PS_W_KNIGHT = 2 * SQUARE_NB,
        PS_B_KNIGHT = 3 * SQUARE_NB,
        PS_W_BISHOP = 4 * SQUARE_NB,
        PS_B_BISHOP = 5 * SQUARE_NB,
        PS_W_ROOK   = 6 * SQUARE_NB,
        PS_B_ROOK   = 7 * SQUARE_NB,
        PS_W_QUEEN  = 8 * SQUARE_NB,
        PS_B_QUEEN  = 9 * SQUARE_NB,
        PS_KING     = 10 * SQUARE_NB,
        PS_NB       = 11 * SQUARE_NB
    };

    static constexpr IndexType PieceSquareIndex[COLOR_NB][PIECE_NB] = {
      // Convention: W - us, B - them
      // Viewed from other side, W and B are reversed
      {PS_NONE, PS_W_PAWN, PS_W_KNIGHT, PS_W_BISHOP, PS_W_ROOK, PS_W_QUEEN, PS_KING, PS_NONE,
       PS_NONE, PS_B_PAWN, PS_B_KNIGHT, PS_B_BISHOP, PS_B_ROOK, PS_B_QUEEN, PS_KING, PS_NONE},
      {PS_NONE, PS_B_PAWN, PS_B_KNIGHT, PS_B_BISHOP, PS_B_ROOK, PS_B_QUEEN, PS_KING, PS_NONE,
       PS_NONE, PS_W_PAWN, PS_W_KNIGHT, PS_W_BISHOP, PS_W_ROOK, PS_W_QUEEN, PS_KING, PS_NONE}};

   public:
    // Feature name
    static constexpr const char* Name = "HalfKAv2_hm(Friend)";

    // Hash value embedded in the evaluation file
    static constexpr std::uint32_t HashValue = 0x7f234cb8u;

    // Number of feature dimensions
    static constexpr IndexType Dimensions =
      static_cast<IndexType>(SQUARE_NB) * static_cast<IndexType>(PS_NB) / 2;

#define B(v) (v * PS_NB)
#define K(f, r) (((f) >= FILE_E ? (7 - (f)) : (f)) * 4 + ((r) > 3 ? (7 - (r)) : (r)))
    // clang-format off
    static constexpr int KingBuckets[COLOR_NB][SQUARE_NB] = {
      { B(K(0,7)), B(K(1,7)), B(K(2,7)), B(K(3,7)), B(K(4,7)), B(K(5,7)), B(K(6,7)), B(K(7,7)),
        B(K(0,6)), B(K(1,6)), B(K(2,6)), B(K(3,6)), B(K(4,6)), B(K(5,6)), B(K(6,6)), B(K(7,6)),
        B(K(0,5)), B(K(1,5)), B(K(2,5)), B(K(3,5)), B(K(4,5)), B(K(5,5)), B(K(6,5)), B(K(7,5)),
        B(K(0,4)), B(K(1,4)), B(K(2,4)), B(K(3,4)), B(K(4,4)), B(K(5,4)), B(K(6,4)), B(K(7,4)),
        B(K(0,3)), B(K(1,3)), B(K(2,3)), B(K(3,3)), B(K(4,3)), B(K(5,3)), B(K(6,3)), B(K(7,3)),
        B(K(0,2)), B(K(1,2)), B(K(2,2)), B(K(3,2)), B(K(4,2)), B(K(5,2)), B(K(6,2)), B(K(7,2)),
        B(K(0,1)), B(K(1,1)), B(K(2,1)), B(K(3,1)), B(K(4,1)), B(K(5,1)), B(K(6,1)), B(K(7,1)),
        B(K(0,0)), B(K(1,0)), B(K(2,0)), B(K(3,0)), B(K(4,0)), B(K(5,0)), B(K(6,0)), B(K(7,0)) },
      { B(K(0,0)), B(K(1,0)), B(K(2,0)), B(K(3,0)), B(K(4,0)), B(K(5,0)), B(K(6,0)), B(K(7,0)),
        B(K(0,1)), B(K(1,1)), B(K(2,1)), B(K(3,1)), B(K(4,1)), B(K(5,1)), B(K(6,1)), B(K(7,1)),
        B(K(0,2)), B(K(1,2)), B(K(2,2)), B(K(3,2)), B(K(4,2)), B(K(5,2)), B(K(6,2)), B(K(7,2)),
        B(K(0,3)), B(K(1,3)), B(K(2,3)), B(K(3,3)), B(K(4,3)), B(K(5,3)), B(K(6,3)), B(K(7,3)),
        B(K(0,4)), B(K(1,4)), B(K(2,4)), B(K(3,4)), B(K(4,4)), B(K(5,4)), B(K(6,4)), B(K(7,4)),
        B(K(0,5)), B(K(1,5)), B(K(2,5)), B(K(3,5)), B(K(4,5)), B(K(5,5)), B(K(6,5)), B(K(7,5)),
        B(K(0,6)), B(K(1,6)), B(K(2,6)), B(K(3,6)), B(K(4,6)), B(K(5,6)), B(K(6,6)), B(K(7,6)),
        B(K(0,7)), B(K(1,7)), B(K(2,7)), B(K(3,7)), B(K(4,7)), B(K(5,7)), B(K(6,7)), B(K(7,7)) }
    };
    // clang-format on
#undef K
#undef B
    // clang-format off
    // Orient a square according to perspective (rotates by 180 for black)
    static constexpr int OrientTBL[COLOR_NB][SQUARE_NB] = {
      { SQ_H1, SQ_H1, SQ_H1, SQ_H1, SQ_A1, SQ_A1, SQ_A1, SQ_A1,
        SQ_H1, SQ_H1, SQ_H1, SQ_H1, SQ_A1, SQ_A1, SQ_A1, SQ_A1,
        SQ_H1, SQ_H1, SQ_H1, SQ_H1, SQ_A1, SQ_A1, SQ_A1, SQ_A1,
        SQ_H1, SQ_H1, SQ_H1, SQ_H1, SQ_A1, SQ_A1, SQ_A1, SQ_A1,
        SQ_H1, SQ_H1, SQ_H1, SQ_H1, SQ_A1, SQ_A1, SQ_A1, SQ_A1,
        SQ_H1, SQ_H1, SQ_H1, SQ_H1, SQ_A1, SQ_A1, SQ_A1, SQ_A1,
        SQ_H1, SQ_H1, SQ_H1, SQ_H1, SQ_A1, SQ_A1, SQ_A1, SQ_A1,
        SQ_H1, SQ_H1, SQ_H1, SQ_H1, SQ_A1, SQ_A1, SQ_A1, SQ_A1 },
      { SQ_H8, SQ_H8, SQ_H8, SQ_H8, SQ_A8, SQ_A8, SQ_A8, SQ_A8,
        SQ_H8, SQ_H8, SQ_H8, SQ_H8, SQ_A8, SQ_A8, SQ_A8, SQ_A8,
        SQ_H8, SQ_H8, SQ_H8, SQ_H8, SQ_A8, SQ_A8, SQ_A8, SQ_A8,
        SQ_H8, SQ_H8, SQ_H8, SQ_H8, SQ_A8, SQ_A8, SQ_A8, SQ_A8,
        SQ_H8, SQ_H8, SQ_H8, SQ_H8, SQ_A8, SQ_A8, SQ_A8, SQ_A8,
        SQ_H8, SQ_H8, SQ_H8, SQ_H8, SQ_A8, SQ_A8, SQ_A8, SQ_A8,
        SQ_H8, SQ_H8, SQ_H8, SQ_H8, SQ_A8, SQ_A8, SQ_A8, SQ_A8,
        SQ_H8, SQ_H8, SQ_H8, SQ_H8, SQ_A8, SQ_A8, SQ_A8, SQ_A8 }
    };
    // clang-format on

    // Maximum number of simultaneously active features.
    static constexpr IndexType MaxActiveDimensions = 32;
    using IndexList                                = ValueList<IndexType, MaxActiveDimensions>;

    // Index of a feature for a given king position and another piece on some square
    template<Color Perspective>
    static IndexType make_index(Square s, Piece pc, Square ksq);

    // Get a list of indices for active features
    template<Color Perspective>
    static void append_active_indices(const Position& pos, IndexList& active);

    // Get a list of indices for recently changed features
    template<Color Perspective>
    static void
    append_changed_indices(Square ksq, const DirtyPiece& dp, IndexList& removed, IndexList& added);

    // Returns the cost of updating one perspective, the most costly one.
    // Assumes no refresh needed.
    static int update_cost(const StateInfo* st);
    static int refresh_cost(const Position& pos);

    // Returns whether the change stored in this StateInfo means
    // that a full accumulator refresh is required.
    static bool requires_refresh(const StateInfo* st, Color perspective);
};

}  // namespace Stockfish::Eval::NNUE::Features

#endif  // #ifndef NNUE_FEATURES_HALF_KA_V2_HM_H_INCLUDED
