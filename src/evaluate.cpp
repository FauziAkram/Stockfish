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
#include <sstream>

#include "nnue/network.h"
#include "nnue/nnue_misc.h"
#include "position.h"
#include "types.h"
#include "uci.h"

namespace Stockfish {
int ee1=500, ee2=513, ee3=517, ee4=499, ee5=32395, ee6=32857, ee7=32793,
    ee8=64, ee9=64, ee10=64, ee11=919, ee12=908, ee13=903, ee14=11, ee15=7,
    ee16=9, ee17=145, ee18=155, ee19=147, ee20=1036 ,ee21=1019 ,ee22=1067,
    ee23=178, ee24=224, ee25=208, ee26=204, ee27=238, ee28=211;
TUNE(SetRange(1, 1501), ee1,ee2,ee3,ee4);
TUNE(SetRange(1, 67001), ee5,ee6,ee7);
TUNE(SetRange(1, 149), ee8,ee9,ee10);
TUNE(SetRange(1, 2001), ee11,ee12,ee13);
TUNE(SetRange(1, 25), ee14,ee15,ee16);
TUNE(SetRange(1, 371), ee17,ee18,ee19);
TUNE(SetRange(1, 2300), ee20,ee21,ee22);
TUNE(SetRange(1, 491), ee23,ee24,ee25);
TUNE(SetRange(1, 521), ee26,ee27,ee28);
// Returns a static, purely materialistic evaluation of the position from
// the point of view of the given color. It can be divided by PawnValue to get
// an approximation of the material advantage on the board in terms of pawns.
int Eval::simple_eval(const Position& pos, Color c) {
    return PawnValue * (pos.count<PAWN>(c) - pos.count<PAWN>(~c))
         + (pos.non_pawn_material(c) - pos.non_pawn_material(~c));
}


// Evaluate is the evaluator for the outer world. It returns a static evaluation
// of the position from the point of view of the side to move.
Value Eval::evaluate(const Eval::NNUE::Networks& networks, const Position& pos, int optimism) {

    assert(!pos.checkers());

    int  simpleEval = simple_eval(pos, pos.side_to_move());
    bool smallNet   = std::abs(simpleEval) > SmallNetThreshold;
    bool psqtOnly   = std::abs(simpleEval) > PsqtOnlyThreshold;
    int  nnueComplexity;
    int  v;

    Value nnue = smallNet ? networks.small.evaluate(pos, true, &nnueComplexity, psqtOnly)
                          : networks.big.evaluate(pos, true, &nnueComplexity, false);

    const auto adjustEval = [&](int optDiv, int nnueDiv, int npmDiv, int pawnCountConstant, int pawnCountMul,
                                int npmConstant, int evalDiv, int shufflingConstant,
                                int shufflingDiv) {
        // Blend optimism and eval with nnue complexity and material imbalance
        optimism += optimism * (nnueComplexity + std::abs(simpleEval - nnue)) / optDiv;
        nnue -= nnue * (nnueComplexity * ee1 / 300) / nnueDiv;

        int npm = pos.non_pawn_material() / npmDiv;
        v       = (nnue * (npm + pawnCountConstant + pawnCountMul * pos.count<PAWN>())
             + optimism * (npmConstant + npm))
          / evalDiv;

        // Damp down the evaluation linearly when shuffling
        int shuffling = pos.rule50_count();
        v             = v * (shufflingConstant - shuffling) / shufflingDiv;
    };

    if (!smallNet)
        adjustEval(ee2, ee5, ee8, ee11, ee14, ee17, ee20, ee23, ee26);
    else if (psqtOnly)
        adjustEval(ee3, ee6, ee9, ee12, ee15, ee18, ee21, ee24, ee27);
    else
        adjustEval(ee4, ee7, ee10, ee13, ee16, ee19, ee22, ee25, ee28);

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

    std::stringstream ss;
    ss << std::showpoint << std::noshowpos << std::fixed << std::setprecision(2);
    ss << '\n' << NNUE::trace(pos, networks) << '\n';

    ss << std::showpoint << std::showpos << std::fixed << std::setprecision(2) << std::setw(15);

    Value v = networks.big.evaluate(pos, false);
    v       = pos.side_to_move() == WHITE ? v : -v;
    ss << "NNUE evaluation        " << 0.01 * UCIEngine::to_cp(v, pos) << " (white side)\n";

    v = evaluate(networks, pos, VALUE_ZERO);
    v = pos.side_to_move() == WHITE ? v : -v;
    ss << "Final evaluation       " << 0.01 * UCIEngine::to_cp(v, pos) << " (white side)";
    ss << " [with scaled NNUE, ...]";
    ss << "\n";

    return ss.str();
}

}  // namespace Stockfish
