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

#include "bitboard.h"
#include "nnue/network.h"
#include "nnue/nnue_misc.h"
#include "position.h"
#include "types.h"
#include "uci.h"
#include "nnue/nnue_accumulator.h"
#include "tune.h"

namespace Stockfish {

namespace Eval {

namespace HCE {

// This namespace holds all the tunable parameters for our
// Hand-Crafted Evaluation (HCE) layer. All are initialized to 0.

// == Piece-specific bonuses ==
int BishopPairBonus = 0;
int KnightOutpostBonus = 0;
int RookOnOpenFile = 0;
int RookOnSemiOpenFile = 0;

// == Pawn structure penalties/bonuses ==
int PassedPawnBonus[RANK_NB] = { 0 };
int IsolatedPawnPenalty = 0;
int DoubledPawnPenalty = 0;
int BackwardPawnPenalty = 0;

// == King safety ==
int KingPawnShieldBonus = 0;
int KingAttackWeight[PIECE_TYPE_NB] = { 0 }; // Weight for each piece type attacking the king zone

// == Mobility ==
int MobilityBonus[PIECE_TYPE_NB] = { 0 }; // Bonus per square of mobility for each piece type

// == Register all parameters for tuning ==
TUNE(SetRange(-200, 200), BishopPairBonus, KnightOutpostBonus);
TUNE(SetRange(-150, 250), RookOnOpenFile, RookOnSemiOpenFile);
TUNE(SetRange(-150, 350), PassedPawnBonus);
TUNE(SetRange(-150, 150), IsolatedPawnPenalty, DoubledPawnPenalty, BackwardPawnPenalty);
TUNE(SetRange(-150, 250), KingPawnShieldBonus);
TUNE(SetRange(-150, 200), KingAttackWeight);
TUNE(SetRange(-120, 120), MobilityBonus);

} // namespace HCE

// hce_evaluate computes the total score from all handcrafted features.
// The score is returned from White's point of view.
Value hce_evaluate(const Position& pos) {
    Value score = 0;

    Bitboard whitePawns = pos.pieces(WHITE, PAWN);
    Bitboard blackPawns = pos.pieces(BLACK, PAWN);
    Bitboard allPawns = whitePawns | blackPawns;

    // 1. Bishop Pair
    if (pos.count<BISHOP>(WHITE) >= 2) score += HCE::BishopPairBonus;
    if (pos.count<BISHOP>(BLACK) >= 2) score -= HCE::BishopPairBonus;

    // 2. Pawns (Passed, Isolated, Doubled, Backward)
    for (Color c : {WHITE, BLACK})
    {
        Bitboard pawns = pos.pieces(c, PAWN);
        Bitboard oppPawns = pos.pieces(~c, PAWN);
        Direction push = pawn_push(c);
        int sign = (c == WHITE ? 1 : -1);

        for (Bitboard b = pawns; b; )
        {
            Square s = pop_lsb(b);
            File f = file_of(s);
            Rank r_real = rank_of(s);
            Rank r_rel = relative_rank(c, s);

            // Passed Pawns: No enemy pawns in front on same or adjacent files
            Bitboard forwardRanks = 0;
            if (c == WHITE) {
                for (Rank r_idx = Rank(r_real + 1); r_idx < RANK_NB; ++r_idx) forwardRanks |= rank_bb(r_idx);
            } else {
                for (Rank r_idx = Rank(r_real - 1); r_idx >= RANK_1; --r_idx) forwardRanks |= rank_bb(r_idx);
            }
            Bitboard passedMask = (file_bb(f) | (shift<WEST>(file_bb(f)) | shift<EAST>(file_bb(f)))) & forwardRanks;
            if (!(passedMask & oppPawns))
            {
                score += sign * HCE::PassedPawnBonus[r_rel];
            }

            // Isolated Pawns
            if (!((shift<WEST>(file_bb(f)) | shift<EAST>(file_bb(f))) & pawns))
            {
                score -= sign * HCE::IsolatedPawnPenalty;
            }

            // Doubled Pawns
            if (popcount(file_bb(f) & pawns) > 1)
            {
                score -= sign * HCE::DoubledPawnPenalty;
            }
            
            // Backward Pawns
            Square stop = s + push;
            if (relative_rank(c, stop) < RANK_8 && !(pawns & attacks_bb<PAWN>(stop, c)))
            {
                if (attacks_bb<PAWN>(stop, ~c) & oppPawns)
                {
                     score -= sign * HCE::BackwardPawnPenalty;
                }
            }
        }
    }

    // 3. Piece Activity (Rooks on open files, Knight outposts)
    for (Color c : {WHITE, BLACK})
    {
        int sign = (c == WHITE ? 1 : -1);

        // Rooks
        for (Bitboard b = pos.pieces(c, ROOK); b; )
        {
            Square s = pop_lsb(b);
            File f = file_of(s);
            if (!(file_bb(f) & allPawns))
                score += sign * HCE::RookOnOpenFile;
            else if (!(file_bb(f) & pos.pieces(c, PAWN)))
                score += sign * HCE::RookOnSemiOpenFile;
        }

        // Knights
        for (Bitboard b = pos.pieces(c, KNIGHT); b; )
        {
            Square s = pop_lsb(b);
            // An outpost is a square supported by a friendly pawn and not attackable by an enemy pawn.
            if ((attacks_bb<PAWN>(s, ~c) & pos.pieces(c, PAWN)) && !(attacks_bb<PAWN>(s, c) & pos.pieces(~c, PAWN)))
            {
                score += sign * HCE::KnightOutpostBonus;
            }
        }
    }
    
    // 4. King Safety (Pawn shield and attacks on king)
    for (Color c : {WHITE, BLACK})
    {
        int sign = (c == WHITE ? 1 : -1);
        Square ksq = pos.square<KING>(c);
        Bitboard kingZone = attacks_bb<KING>(ksq) | attacks_bb<KING>(ksq + pawn_push(c));
        
        // Pawn Shield
        score += sign * popcount(kingZone & pos.pieces(c, PAWN)) * HCE::KingPawnShieldBonus;
        
        // Attacks on King
        for (PieceType pt = KNIGHT; pt <= QUEEN; ++pt)
        {
            Bitboard attackers = pos.pieces(~c, pt);
            while (attackers)
            {
                Square s = pop_lsb(attackers);
                score -= sign * popcount(attacks_bb(pt, s, pos.pieces()) & kingZone) * HCE::KingAttackWeight[pt];
            }
        }
    }

    // 5. Mobility
    for (Color c : {WHITE, BLACK})
    {
        int sign = (c == WHITE ? 1 : -1);
        for (PieceType pt = KNIGHT; pt <= QUEEN; ++pt)
        {
            Bitboard pieces = pos.pieces(c, pt);
            while (pieces)
            {
                Square s = pop_lsb(pieces);
                score += sign * popcount(attacks_bb(pt, s, pos.pieces()) & ~pos.pieces(c)) * HCE::MobilityBonus[pt];
            }
        }
    }

    return score;
}

} // namespace Eval

// Returns a static, purely materialistic evaluation of the position from
// the point of view of the side to move. It can be divided by PawnValue to get
// an approximation of the material advantage on the board in terms of pawns.
int Eval::simple_eval(const Position& pos) {
    Color c = pos.side_to_move();
    return PawnValue * (pos.count<PAWN>(c) - pos.count<PAWN>(~c))
         + (pos.non_pawn_material(c) - pos.non_pawn_material(~c));
}

bool Eval::use_smallnet(const Position& pos) { return std::abs(simple_eval(pos)) > 962; }

// Evaluate is the evaluator for the outer world. It returns a static evaluation
// of the position from the point of view of the side to move.
Value Eval::evaluate(const Eval::NNUE::Networks&    networks,
                     const Position&                pos,
                     Eval::NNUE::AccumulatorStack&  accumulators,
                     Eval::NNUE::AccumulatorCaches& caches,
                     int                            optimism) {

    assert(!pos.checkers());

    bool smallNet           = use_smallnet(pos);
    auto [psqt, positional] = smallNet ? networks.small.evaluate(pos, accumulators, &caches.small)
                                       : networks.big.evaluate(pos, accumulators, &caches.big);

    Value nnue = (125 * psqt + 131 * positional) / 128;

    // HCE Layer: Add handcrafted bonuses/penalties on top of NNUE evaluation
    Value hce_bonus_w_pov = hce_evaluate(pos);
    nnue += (pos.side_to_move() == WHITE) ? hce_bonus_w_pov : -hce_bonus_w_pov;

    // Re-evaluate the position when higher eval accuracy is worth the time spent
    if (smallNet && (std::abs(nnue) < 236))
    {
        std::tie(psqt, positional) = networks.big.evaluate(pos, accumulators, &caches.big);
        nnue                       = (125 * psqt + 131 * positional) / 128;
        nnue += (pos.side_to_move() == WHITE) ? hce_bonus_w_pov : -hce_bonus_w_pov; // Re-apply HCE
        smallNet                   = false;
    }

    // Blend optimism and eval with nnue complexity
    int nnueComplexity = std::abs(psqt - positional);
    optimism += optimism * nnueComplexity / 468;
    nnue -= nnue * nnueComplexity / 18000;

    int material = 535 * pos.count<PAWN>() + pos.non_pawn_material();
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

    Eval::NNUE::AccumulatorStack accumulators;
    auto                         caches = std::make_unique<Eval::NNUE::AccumulatorCaches>(networks);

    std::stringstream ss;
    ss << std::showpoint << std::noshowpos << std::fixed << std::setprecision(2);
    ss << '\n' << NNUE::trace(pos, networks, *caches) << '\n';

    ss << std::showpoint << std::showpos << std::fixed << std::setprecision(2) << std::setw(15);

    auto [psqt, positional] = networks.big.evaluate(pos, accumulators, &caches->big);
    Value v                 = psqt + positional;
    v                       = pos.side_to_move() == WHITE ? v : -v;
    ss << "NNUE evaluation        " << 0.01 * UCIEngine::to_cp(v, pos) << " (white side)\n";

    Value hce_bonus_w_pov = hce_evaluate(pos);
    ss << "HCE bonus              " << 0.01 * UCIEngine::to_cp(hce_bonus_w_pov, pos) << " (white side)\n";

    v = evaluate(networks, pos, accumulators, *caches, 0);
    v = pos.side_to_move() == WHITE ? v : -v;
    ss << "Final evaluation       " << 0.01 * UCIEngine::to_cp(v, pos) << " (white side)";
    ss << " [with HCE, scaled NNUE, ...]";
    ss << "\n";

    return ss.str();
}

}  // namespace Stockfish
