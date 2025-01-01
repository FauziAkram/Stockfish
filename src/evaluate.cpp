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

int AA1=1048, AA2=1048,AA3=1048, AA4=1048,AA5=1048, AA6=1048,AA7=1048, AA8=1048,AA9=1048;
int BB1=1048, BB2=1048,BB3=1048, BB4=1048,BB5=1048, BB6=1048,BB7=1048, BB8=1048,BB9=1048;
int CC1=1048, CC2=1048,CC3=1048, CC4=1048,CC5=1048, CC6=1048,CC7=1048, CC8=1048,CC9=1048;
int DD1=1048, DD2=1048,DD3=1048, DD4=1048,DD5=1048, DD6=1048,DD7=1048, DD8=1048,DD9=1048;
int EE1=1048, EE2=1048,EE3=1048, EE4=1048,EE5=1048, EE6=1048,EE7=1048, EE8=1048,EE9=1048;
int FF1=1048, FF2=1048,FF3=1048, FF4=1048,FF5=1048, FF6=1048,FF7=1048, FF8=1048,FF9=1048;
int GG1=1048, GG2=1048,GG3=1048, GG4=1048,GG5=1048, GG6=1048,GG7=1048, GG8=1048,GG9=1048;
int HH1=1048, HH2=1048,HH3=1048, HH4=1048,HH5=1048, HH6=1048,HH7=1048, HH8=1048,HH9=1048;
int II1=1048, II2=1048,II3=1048, II4=1048,II5=1048, II6=1048,II7=1048, II8=1048,II9=1048;

TUNE(AA1, AA2,AA3, AA4,AA5, AA6,AA7, AA8,AA9);
TUNE(BB1, BB2,BB3, BB4,BB5, BB6,BB7, BB8,BB9);
TUNE(CC1, CC2,CC3, CC4,CC5, CC6,CC7, CC8,CC9);
TUNE(DD1, DD2,DD3, DD4,DD5, DD6,DD7, DD8,DD9);
TUNE(EE1, EE2,EE3, EE4,EE5, EE6,EE7, EE8,EE9);
TUNE(FF1, FF2,FF3, FF4,FF5, FF6,FF7, FF8,FF9);
TUNE(GG1, GG2,GG3, GG4,GG5, GG6,GG7, GG8,GG9);
TUNE(HH1, HH2,HH3, HH4,HH5, HH6,HH7, HH8,HH9);
TUNE(II1, II2,II3, II4,II5, II6,II7, II8,II9);


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

    Value scaled_psqt = 1000 * psqt;
  
  Value scaled_positional = ((psqt < -3500)?  ((positional < -3500)? AA1:
                                     (positional < -2500)? AA2:
                                     (positional < -1500)? AA3:
                                     (positional < -500)? AA4:
                                     (positional < 500)? AA5:
                                     (positional < 1500)? AA6:
                                     (positional < 2500)? AA7:
                                     (positional < 3500)? AA8:AA9):
                     (psqt < -2500)?  ((positional < -3500)? BB1:
                                     (positional < -2500)? BB2:
                                     (positional < -1500)? BB3:
                                     (positional < -500)? BB4:
                                     (positional < 500)? BB5:
                                     (positional < 1500)? BB6:
                                     (positional < 2500)? BB7:
                                     (positional < 3500)? BB8:BB9):
                     (psqt < -1500)?  ((positional < -3500)? CC1:
                                     (positional < -2500)? CC2:
                                     (positional < -1500)? CC3:
                                     (positional < -500)? CC4:
                                     (positional < 500)? CC5:
                                     (positional < 1500)? CC6:
                                     (positional < 2500)? CC7:
                                     (positional < 3500)? CC8:CC9):
                     (psqt < -500)?  ((positional < -3500)? DD1:
                                     (positional < -2500)? DD2:
                                     (positional < -1500)? DD3:
                                     (positional < -500)? DD4:
                                     (positional < 500)? DD5:
                                     (positional < 1500)? DD6:
                                     (positional < 2500)? DD7:
                                     (positional < 3500)? DD8:DD9):
                     (psqt < 500)?  ((positional < -3500)? EE1:
                                     (positional < -2500)? EE2:
                                     (positional < -1500)? EE3:
                                     (positional < -500)? EE4:
                                     (positional < 500)? EE5:
                                     (positional < 1500)? EE6:
                                     (positional < 2500)? EE7:
                                     (positional < 3500)? EE8:EE9):
                    (psqt < 1500)?  ((positional < -3500)? FF1:
                                     (positional < -2500)? FF2:
                                     (positional < -1500)? FF3:
                                     (positional < -500)? FF4:
                                     (positional < 500)? FF5:
                                     (positional < 1500)? FF6:
                                     (positional < 2500)? FF7:
                                     (positional < 3500)? FF8:FF9):
                    (psqt < 2500)?  ((positional < -3500)? GG1:
                                     (positional < -2500)? GG2:
                                     (positional < -1500)? GG3:
                                     (positional < -500)? GG4:
                                     (positional < 500)? GG5:
                                     (positional < 1500)? GG6:
                                     (positional < 2500)? GG7:
                                     (positional < 3500)? GG8:GG9):
                     (psqt < 3500)?  ((positional < -3500)? HH1:
                                     (positional < -2500)? HH2:
                                     (positional < -1500)? HH3:
                                     (positional < -500)? HH4:
                                     (positional < 500)? HH5:
                                     (positional < 1500)? HH6:
                                     (positional < 2500)? HH7:
                                     (positional < 3500)? HH8:HH9):
                                    ((positional < -3500)? II1:
                                     (positional < -2500)? II2:
                                     (positional < -1500)? II3:
                                     (positional < -500)? II4:
                                     (positional < 500)? II5:
                                     (positional < 1500)? II6:
                                     (positional < 2500)? II7:
                                     (positional < 3500)? II8:II9)) * positional;

    Value nnue = (scaled_psqt + scaled_positional) / 1024;

    // Re-evaluate the position when higher eval accuracy is worth the time spent
    if (smallNet && (std::abs(nnue) < 236))
    {
        std::tie(psqt, positional) = networks.big.evaluate(pos, &caches.big);
        nnue                       = (1000 * psqt + 1048 * positional) / 1024;
        smallNet                   = false;
    }

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
