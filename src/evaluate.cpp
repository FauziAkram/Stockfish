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
int aa1=1000, aa2=1000,aa3=1000, aa4=1000,aa5=1000, aa6=1000,aa7=1000, aa8=1000,aa9=1000, aa10=1000;
int bb1=1000, bb2=1000,bb3=1000, bb4=1000,bb5=1000, bb6=1000,bb7=1000, bb8=1000,bb9=1000, bb10=1000;
int cc1=1000, cc2=1000,cc3=1000, cc4=1000,cc5=1000, cc6=1000,cc7=1000, cc8=1000,cc9=1000, cc10=1000;
int dd1=1000, dd2=1000,dd3=1000, dd4=1000,dd5=1000, dd6=1000,dd7=1000, dd8=1000,dd9=1000, dd10=1000;
int ee1=1000, ee2=1000,ee3=1000, ee4=1000,ee5=1000, ee6=1000,ee7=1000, ee8=1000,ee9=1000, ee10=1000;
int ff1=1000, ff2=1000,ff3=1000, ff4=1000,ff5=1000, ff6=1000,ff7=1000, ff8=1000,ff9=1000, ff10=1000;
int gg1=1000, gg2=1000,gg3=1000, gg4=1000,gg5=1000, gg6=1000,gg7=1000, gg8=1000,gg9=1000, gg10=1000;
int hh1=1000, hh2=1000,hh3=1000, hh4=1000,hh5=1000, hh6=1000,hh7=1000, hh8=1000,hh9=1000, hh10=1000;
int ii1=1000, ii2=1000,ii3=1000, ii4=1000,ii5=1000, ii6=1000,ii7=1000, ii8=1000,ii9=1000, ii10=1000;

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

    Value nnue = (  (psqt < -3500)?  ((positional < -3500)? aa1:
                                     (positional < -2500)? aa2:
                                     (positional < -1500)? aa3:
                                     (positional < -500)? aa4:
                                     (positional < 500)? aa5:
                                     (positional < 1500)? aa6:
                                     (positional < 2500)? aa7:
                                     (positional < 3500)? aa8:aa9):
                     (psqt < -2500)?  ((positional < -3500)? bb1:
                                     (positional < -2500)? bb2:
                                     (positional < -1500)? bb3:
                                     (positional < -500)? bb4:
                                     (positional < 500)? bb5:
                                     (positional < 1500)? bb6:
                                     (positional < 2500)? bb7:
                                     (positional < 3500)? bb8:bb9):
                     (psqt < -1500)?  ((positional < -3500)? cc1:
                                     (positional < -2500)? cc2:
                                     (positional < -1500)? cc3:
                                     (positional < -500)? cc4:
                                     (positional < 500)? cc5:
                                     (positional < 1500)? cc6:
                                     (positional < 2500)? cc7:
                                     (positional < 3500)? cc8:cc9):
                     (psqt < -500)?  ((positional < -3500)? dd1:
                                     (positional < -2500)? dd2:
                                     (positional < -1500)? dd3:
                                     (positional < -500)? dd4:
                                     (positional < 500)? dd5:
                                     (positional < 1500)? dd6:
                                     (positional < 2500)? dd7:
                                     (positional < 3500)? dd8:dd9):
                     (psqt < 500)?  ((positional < -3500)? ee1:
                                     (positional < -2500)? ee2:
                                     (positional < -1500)? ee3:
                                     (positional < -500)? ee4:
                                     (positional < 500)? ee5:
                                     (positional < 1500)? ee6:
                                     (positional < 2500)? ee7:
                                     (positional < 3500)? ee8:ee9):
                    (psqt < 1500)?  ((positional < -3500)? ff1:
                                     (positional < -2500)? ff2:
                                     (positional < -1500)? ff3:
                                     (positional < -500)? ff4:
                                     (positional < 500)? ff5:
                                     (positional < 1500)? ff6:
                                     (positional < 2500)? ff7:
                                     (positional < 3500)? ff8:ff9):
                    (psqt < 2500)?  ((positional < -3500)? gg1:
                                     (positional < -2500)? gg2:
                                     (positional < -1500)? gg3:
                                     (positional < -500)? gg4:
                                     (positional < 500)? gg5:
                                     (positional < 1500)? gg6:
                                     (positional < 2500)? gg7:
                                     (positional < 3500)? gg8:gg9):
                     (psqt < 3500)?  ((positional < -3500)? hh1:
                                     (positional < -2500)? hh2:
                                     (positional < -1500)? hh3:
                                     (positional < -500)? hh4:
                                     (positional < 500)? hh5:
                                     (positional < 1500)? hh6:
                                     (positional < 2500)? hh7:
                                     (positional < 3500)? hh8:hh9):
                                    ((positional < -3500)? ii1:
                                     (positional < -2500)? ii2:
                                     (positional < -1500)? ii3:
                                     (positional < -500)? ii4:
                                     (positional < 500)? ii5:
                                     (positional < 1500)? ii6:
                                     (positional < 2500)? ii7:
                                     (positional < 3500)? ii8:ii9)
    * psqt + 1048 * positional) / 1024;

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
