/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2024 The Stockfish developers (see AUTHORS file)

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
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <tuple>

#include "nnue/network.h"
#include "nnue/nnue_misc.h"
#include "position.h"
#include "types.h"
#include "uci.h"
#include "nnue/nnue_accumulator.h"

namespace Stockfish {

// Returns a static, purely materialistic evaluation of the position from
// the point of view of the given color. It can be divided by PawnValue to get
// an approximation of the material advantage on the board in terms of pawns.
int Eval::simple_eval(const Position& pos, Color c) {
    return PawnValue * (pos.count<PAWN>(c) - pos.count<PAWN>(~c))
         + (pos.non_pawn_material(c) - pos.non_pawn_material(~c));
}

bool Eval::use_smallnet(const Position& pos) {
    int simpleEval = simple_eval(pos, pos.side_to_move());
    return std::abs(simpleEval) > 962;
}

// Evaluate is the evaluator for the outer world. It returns a static evaluation
// of the position from the point of view of the side to move.
Value Eval::evaluate(const Eval::NNUE::Networks&    networks,
                     const Position&                pos,
                     Eval::NNUE::AccumulatorCaches& caches,
                     int                            optimism) {

    assert(!pos.checkers());

    bool smallNet           = use_smallnet(pos);
    auto [psqt, positional] = smallNet ? networks.small.evaluate(pos, &caches.small)
                                       : networks.big.evaluate(pos, &caches.big);

    Value nnue = (125 * psqt + 131 * positional) / 128;

    // Re-evaluate the position when higher eval accuracy is worth the time spent
    if (smallNet && (std::abs(nnue) < 236))
    {
        std::tie(psqt, positional) = networks.big.evaluate(pos, &caches.big);
        nnue                       = (125 * psqt + 131 * positional) / 128;
        smallNet                   = false;
    }

    // 1. Threatened Pawns and Safe Mobility
    Bitboard threatenedByPawn = pos.attacks_by<PAWN>(~pos.side_to_move());
    Bitboard safeSquares = ~pos.pieces() & ~threatenedByPawn;

    Bitboard pawnMobility = 0;
    if (pos.side_to_move() == WHITE)
        pawnMobility = (pos.pieces(WHITE, PAWN) << 8) & safeSquares;
    else
        pawnMobility = (pos.pieces(BLACK, PAWN) >> 8) & safeSquares;

    int mobilityBonus = popcount(pawnMobility);
  
    // 2. Isolated Pawns
    int isolatedPawnPenalty = 0;
    Bitboard ourPawns = pos.pieces(pos.side_to_move(), PAWN);
    // Iterate using pop_lsb
    Bitboard pawnsCopy = ourPawns; // Work on a copy to avoid modifying the original
    while (pawnsCopy) {
        Square s = pop_lsb(pawnsCopy); // Get and remove the least significant bit
        File f = file_of(s);
        Bitboard adjacentFiles = (f != FILE_A ? file_bb(File(f - 1)) : 0) | (f != FILE_H ? file_bb(File(f + 1)) : 0);
        if (!(adjacentFiles & ourPawns))
            isolatedPawnPenalty--;
    }

    // 3. Doubled Pawns
    int doubledPawnPenalty = 0;
    for (File f = FILE_A; f <= FILE_H; ++f) {
        Bitboard filePawns = file_bb(f) & ourPawns;
        if (more_than_one(filePawns))
            doubledPawnPenalty -= popcount(filePawns) - 1;
    }

    // 4. Passed Pawns
    int passedPawnBonus = 0;
    Bitboard enemyPawns = pos.pieces(~pos.side_to_move(), PAWN);

    // Iterate using pop_lsb
    pawnsCopy = ourPawns; // Reset copy for reuse
    while (pawnsCopy) {
        Square s = pop_lsb(pawnsCopy);
        File f = file_of(s);
        Rank r = rank_of(s);

        Bitboard path = 0;
        if (pos.side_to_move() == WHITE) {
            for (Rank ahead = Rank(r + 1); ahead <= RANK_8; ++ahead)
                path |= square_bb(make_square(f, ahead));
        } else {
            for (Rank ahead = Rank(r - 1); ahead >= RANK_1; --ahead)
                path |= square_bb(make_square(f, ahead));
        }
        if (f != FILE_A)
            path |= (path << 1);

        if (f != FILE_H)
            path |= (path >> 1);

        if (!(path & enemyPawns))
            passedPawnBonus += (pos.side_to_move() == WHITE ? r : 7 - r); // Based on rank
    }

    // Combine the modifications into the 'positional' component.
    positional += (mobilityBonus * 0) + (isolatedPawnPenalty  * 0) + (doubledPawnPenalty * 0) + (passedPawnBonus * 0);

    nnue = (125 * psqt + 131 * positional) / 128; // Recalculate NNUE with adjusted positional

    // Blend optimism and eval with nnue complexity
    int nnueComplexity = std::abs(psqt - positional);
    optimism += optimism * nnueComplexity / 468;
    nnue -= nnue * nnueComplexity / (smallNet ? 20233 : 17879);

    int material = (smallNet ? 553 : 532) * pos.count<PAWN>() + pos.non_pawn_material();
    int v        = (nnue * (77777 + material) + optimism * (7777 + material)) / 77777;

    // Damp down the evaluation linearly when shuffling
    v -= v * pos.rule50_count() / 212;

    // Guarantee evaluation does not hit the tablebase range
    v = std::clamp(v, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);

    return v;
}

// Like evaluate(), but instead of returning a value, it returns
// a string (suitable for outputting to stdout) that contains the detailed
// descriptions and values of each evaluation term. Useful for debugging.
// Trace scores are from white's point of view
std::string Eval::trace(Position& pos, const Eval::NNUE::Networks& networks) {

    if (pos.checkers())
        return "Final evaluation: none (in check)";

    auto caches = std::make_unique<Eval::NNUE::AccumulatorCaches>(networks);

    std::stringstream ss;
    ss << std::showpoint << std::noshowpos << std::fixed << std::setprecision(2);
    ss << '\n' << NNUE::trace(pos, networks, *caches) << '\n';

    ss << std::showpoint << std::showpos << std::fixed << std::setprecision(2) << std::setw(15);

    auto [psqt, positional] = networks.big.evaluate(pos, &caches->big);
    Value v                 = psqt + positional;
    v                       = pos.side_to_move() == WHITE ? v : -v;
    ss << "NNUE evaluation        " << 0.01 * UCIEngine::to_cp(v, pos) << " (white side)\n";

    v = evaluate(networks, pos, *caches, VALUE_ZERO);
    v = pos.side_to_move() == WHITE ? v : -v;
    ss << "Final evaluation       " << 0.01 * UCIEngine::to_cp(v, pos) << " (white side)";
    ss << " [with scaled NNUE, ...]";
    ss << "\n";

    return ss.str();
}

}  // namespace Stockfish
