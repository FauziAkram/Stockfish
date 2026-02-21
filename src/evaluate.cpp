/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2026 The Stockfish developers (see AUTHORS file)

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

#include "evaluate.h"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <sstream>

#include "bitboard.h"
#include "position.h"
#include "types.h"
#include "uci.h"

namespace Stockfish {

namespace {

constexpr int PieceValue[PIECE_TYPE_NB] = {0, 100, 320, 330, 500, 900, 0};

int center_bonus(Square s) {
    int f = std::abs(file_of(s) - FILE_D);
    int r = std::abs(rank_of(s) - RANK_4);
    return std::max(0, 6 - (f + r));
}

int piece_square_bonus(PieceType pt, Square s, Color c) {
    int relativeRank = c == WHITE ? int(rank_of(s)) : int(RANK_8 - rank_of(s));

    switch (pt)
    {
    case PAWN: return relativeRank * 6 + center_bonus(s);
    case KNIGHT: return center_bonus(s) * 5;
    case BISHOP: return center_bonus(s) * 3;
    case ROOK: return relativeRank * 2;
    case QUEEN: return center_bonus(s) * 2;
    case KING:
        return (relativeRank <= 1 ? 15 : -5 * (relativeRank - 1))
             + (file_of(s) == FILE_C || file_of(s) == FILE_G ? 10 : 0);
    default: return 0;
    }
}

int hce_side_score(const Position& pos, Color c) {
    int score = 0;

    for (PieceType pt : {PAWN, KNIGHT, BISHOP, ROOK, QUEEN, KING})
    {
        score += PieceValue[pt] * popcount(pos.pieces(c, pt));

        Bitboard b = pos.pieces(c, pt);
        while (b)
        {
            Square s = pop_lsb(b);
            score += piece_square_bonus(pt, s, c);
        }
    }

    return score;
}

}  // namespace

int Eval::simple_eval(const Position& pos) {
    Color c = pos.side_to_move();
    return PawnValue * (pos.count<PAWN>(c) - pos.count<PAWN>(~c)) + pos.non_pawn_material(c)
         - pos.non_pawn_material(~c);
}

Value Eval::evaluate(const Position& pos, int optimism) {

    assert(!pos.checkers());

    int hce = hce_side_score(pos, WHITE) - hce_side_score(pos, BLACK);
    int v   = pos.side_to_move() == WHITE ? hce : -hce;

    v += optimism / 8;

    int material = 534 * pos.count<PAWN>() + pos.non_pawn_material();
    v            = v * (63270 + material) / 63270;

    v -= v * pos.rule50_count() / 180;

    v = std::clamp(v, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);

    return Value(v);
}

std::string Eval::trace(Position& pos) {

    if (pos.checkers())
        return "Final evaluation: none (in check)";

    std::stringstream ss;
    ss << std::showpoint << std::showpos << std::fixed << std::setprecision(2) << std::setw(15);

    Value hce = Value(hce_side_score(pos, WHITE) - hce_side_score(pos, BLACK));
    ss << "HCE evaluation         " << 0.01 * UCIEngine::to_cp(hce, pos) << " (white side)\n";

    Value v = evaluate(pos, VALUE_ZERO);
    v       = pos.side_to_move() == WHITE ? v : -v;
    ss << "Final evaluation       " << 0.01 * UCIEngine::to_cp(v, pos) << " (white side)\n";

    return ss.str();
}

}  // namespace Stockfish
