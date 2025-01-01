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
int aa1=1000, aa2=1000,aa3=1000, aa4=1000,aa5=1000, aa6=1000,aa7=1000, aa8=1000,aa9=1000;
int bb1=1000, bb2=1000,bb3=1000, bb4=1000,bb5=1000, bb6=1000,bb7=1000, bb8=1000,bb9=1000;
int cc1=1000, cc2=1000,cc3=1000, cc4=1000,cc5=1000, cc6=1000,cc7=1000, cc8=1000,cc9=1000;
int dd1=1000, dd2=1000,dd3=1000, dd4=1000,dd5=1000, dd6=1000,dd7=1000, dd8=1000,dd9=1000;
int ee1=1000, ee2=1000,ee3=1000, ee4=1000,ee5=1000, ee6=1000,ee7=1000, ee8=1000,ee9=1000;
int ff1=1000, ff2=1000,ff3=1000, ff4=1000,ff5=1000, ff6=1000,ff7=1000, ff8=1000,ff9=1000;
int gg1=1000, gg2=1000,gg3=1000, gg4=1000,gg5=1000, gg6=1000,gg7=1000, gg8=1000,gg9=1000;
int hh1=1000, hh2=1000,hh3=1000, hh4=1000,hh5=1000, hh6=1000,hh7=1000, hh8=1000,hh9=1000;
int ii1=1000, ii2=1000,ii3=1000, ii4=1000,ii5=1000, ii6=1000,ii7=1000, ii8=1000,ii9=1000;
int AA1=1048, AA2=1048,AA3=1048, AA4=1048,AA5=1048, AA6=1048,AA7=1048, AA8=1048,AA9=1048;
int BB1=1048, BB2=1048,BB3=1048, BB4=1048,BB5=1048, BB6=1048,BB7=1048, BB8=1048,BB9=1048;
int CC1=1048, CC2=1048,CC3=1048, CC4=1048,CC5=1048, CC6=1048,CC7=1048, CC8=1048,CC9=1048;
int DD1=1048, DD2=1048,DD3=1048, DD4=1048,DD5=1048, DD6=1048,DD7=1048, DD8=1048,DD9=1048;
int EE1=1048, EE2=1048,EE3=1048, EE4=1048,EE5=1048, EE6=1048,EE7=1048, EE8=1048,EE9=1048;
int FF1=1048, FF2=1048,FF3=1048, FF4=1048,FF5=1048, FF6=1048,FF7=1048, FF8=1048,FF9=1048;
int GG1=1048, GG2=1048,GG3=1048, GG4=1048,GG5=1048, GG6=1048,GG7=1048, GG8=1048,GG9=1048;
int HH1=1048, HH2=1048,HH3=1048, HH4=1048,HH5=1048, HH6=1048,HH7=1048, HH8=1048,HH9=1048;
int II1=1048, II2=1048,II3=1048, II4=1048,II5=1048, II6=1048,II7=1048, II8=1048,II9=1048;
TUNE(aa1, aa2,aa3, aa4,aa5, aa6,aa7, aa8,aa9);
TUNE(bb1, bb2,bb3, bb4,bb5, bb6,bb7, bb8,bb9);
TUNE(cc1, cc2,cc3, cc4,cc5, cc6,cc7, cc8,cc9);
TUNE(dd1, dd2,dd3, dd4,dd5, dd6,dd7, dd8,dd9);
TUNE(ee1, ee2,ee3, ee4,ee5, ee6,ee7, ee8,ee9);
TUNE(ff1, ff2,ff3, ff4,ff5, ff6,ff7, ff8,ff9);
TUNE(gg1, gg2,gg3, gg4,gg5, gg6,gg7, gg8,gg9);
TUNE(hh1, hh2,hh3, hh4,hh5, hh6,hh7, hh8,hh9);
TUNE(ii1, ii2,ii3, ii4,ii5, ii6,ii7, ii8,ii9);
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

    Value nnue = ((  (((psqt < -3500)?  ((positional < -3500)? aa1:
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
                                     (positional < 3500)? ii8:ii9)) * psqt) +
  (((psqt < -3500)?  ((positional < -3500)? AA1:
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
                                     (positional < 3500)? II8:II9)) * positional)) / 1024);

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
