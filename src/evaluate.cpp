/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2026 The Stockfish developers (see AUTHORS file)
*/

#include "evaluate.h"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <sstream>

#include "bitboard.h"
#include "endgame.h"
#include "material.h"
#include "pawns.h"
#include "position.h"
#include "psqt.h"
#include "uci.h"

namespace Stockfish::Eval {

namespace {

int psqt_score(const Position& pos) {
    int score = 0;

    for (Piece pc : {W_PAWN, W_KNIGHT, W_BISHOP, W_ROOK, W_QUEEN, W_KING,
                     B_PAWN, B_KNIGHT, B_BISHOP, B_ROOK, B_QUEEN, B_KING})
    {
        Bitboard b = pos.pieces(type_of(pc));
        b &= pos.pieces(color_of(pc));
        while (b)
            score += PSQT::psq[pc][pop_lsb(b)];
    }

    return score;
}

int mobility_score(const Position& pos, Color c) {
    int ownOcc = popcount(pos.pieces(c));
    int m      = 0;

    for (PieceType pt : {KNIGHT, BISHOP, ROOK, QUEEN})
    {
        Bitboard b = pos.pieces(c, pt);
        while (b)
        {
            Square s = pop_lsb(b);
            Bitboard a = pt == KNIGHT ? attacks_bb<KNIGHT>(s)
                       : pt == BISHOP ? attacks_bb<BISHOP>(s, pos.pieces())
                       : pt == ROOK   ? attacks_bb<ROOK>(s, pos.pieces())
                                      : attacks_bb<QUEEN>(s, pos.pieces());
            m += popcount(a) - ownOcc / 16;
        }
    }

    return m;
}

int passed_score(const Position& pos, const Pawns::Entry* pe, Color c) {
    int bonus = 0;
    Bitboard b = pe->passed_pawns(c);
    while (b)
    {
        Square s = pop_lsb(b);
        int r = int(relative_rank(c, s));
        bonus += r * r * 8;
    }
    return bonus;
}

}  // namespace

Value evaluate(const Position& pos) {

    assert(!pos.checkers());

    Material::Entry* me = Material::probe(pos);
    if (me->specialized_eval_exists())
        return me->evaluate(pos);

    Pawns::Entry* pe = Pawns::probe(pos);

    int score = 0;
    score += psqt_score(pos);
    score += me->imbalance();
    score += pe->pawn_score(WHITE) - pe->pawn_score(BLACK);
    score += pe->king_safety<WHITE>(pos) - pe->king_safety<BLACK>(pos);
    score += passed_score(pos, pe, WHITE) - passed_score(pos, pe, BLACK);
    score += 2 * (mobility_score(pos, WHITE) - mobility_score(pos, BLACK));

    Color strongSide = score >= 0 ? WHITE : BLACK;
    int   sf         = me->scale_factor(pos, strongSide);
    score            = score * sf / SCALE_FACTOR_NORMAL;

    score = (score / 16) * 16;
    score = pos.side_to_move() == WHITE ? score : -score;
    score = score * (200 - pos.rule50_count()) / 214;
    score = std::clamp(score, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);

    return Value(score);
}

std::string trace(Position& pos) {

    if (pos.checkers())
        return "Final evaluation: none (in check)";

    Material::Entry* me = Material::probe(pos);
    Pawns::Entry*    pe = Pawns::probe(pos);

    int psq    = psqt_score(pos);
    int imb    = me->imbalance();
    int pawns  = pe->pawn_score(WHITE) - pe->pawn_score(BLACK);
    int king   = pe->king_safety<WHITE>(pos) - pe->king_safety<BLACK>(pos);
    int passed = passed_score(pos, pe, WHITE) - passed_score(pos, pe, BLACK);
    int mob    = 2 * (mobility_score(pos, WHITE) - mobility_score(pos, BLACK));

    Value v = evaluate(pos);
    Value whiteV = pos.side_to_move() == WHITE ? v : -v;

    std::stringstream ss;
    ss << std::fixed << std::setprecision(2);
    ss << "PSQT                 " << 0.01 * UCIEngine::to_cp(Value(psq), pos) << "\n";
    ss << "Material imbalance   " << 0.01 * UCIEngine::to_cp(Value(imb), pos) << "\n";
    ss << "Pawns                " << 0.01 * UCIEngine::to_cp(Value(pawns), pos) << "\n";
    ss << "King safety          " << 0.01 * UCIEngine::to_cp(Value(king), pos) << "\n";
    ss << "Passed pawns         " << 0.01 * UCIEngine::to_cp(Value(passed), pos) << "\n";
    ss << "Mobility             " << 0.01 * UCIEngine::to_cp(Value(mob), pos) << "\n";
    ss << "Final evaluation     " << 0.01 * UCIEngine::to_cp(whiteV, pos) << " (white side)\n";

    return ss.str();
}

}  // namespace Stockfish::Eval
