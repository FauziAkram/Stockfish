/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
*/

#include "material.h"

#include <algorithm>
#include <cassert>
#include <cstring>

namespace Stockfish {

namespace {
#define S(mg, eg) make_score(mg, eg)

constexpr EvalScore QuadraticOurs[][PIECE_TYPE_NB] = {
  {S(1419, 1455)},
  {S(101, 28), S(37, 39)},
  {S(57, 64), S(249, 187), S(-49, -62)},
  {S(0, 0), S(118, 137), S(10, 27), S(0, 0)},
  {S(-63, -68), S(-5, 3), S(100, 81), S(132, 118), S(-246, -244)},
  {S(-210, -211), S(37, 14), S(147, 141), S(161, 105), S(-158, -174), S(-9, -31)}
};

constexpr EvalScore QuadraticTheirs[][PIECE_TYPE_NB] = {
  {},
  {S(33, 30)},
  {S(46, 18), S(106, 84)},
  {S(75, 35), S(59, 44), S(60, 15)},
  {S(26, 35), S(6, 22), S(38, 39), S(-12, -2)},
  {S(97, 93), S(100, 163), S(-58, -91), S(112, 192), S(276, 225)}
};

#undef S

Endgame<KXK>     EvaluateKXK[] = {Endgame<KXK>(WHITE), Endgame<KXK>(BLACK)};
Endgame<KBPsK>   ScaleKBPsK[]  = {Endgame<KBPsK>(WHITE), Endgame<KBPsK>(BLACK)};
Endgame<KQKRPs>  ScaleKQKRPs[] = {Endgame<KQKRPs>(WHITE), Endgame<KQKRPs>(BLACK)};
Endgame<KPsK>    ScaleKPsK[]   = {Endgame<KPsK>(WHITE), Endgame<KPsK>(BLACK)};
Endgame<KPKP>    ScaleKPKP[]   = {Endgame<KPKP>(WHITE), Endgame<KPKP>(BLACK)};

bool is_KXK(const Position& pos, Color us) {
    return !more_than_one(pos.pieces(~us)) && pos.non_pawn_material(us) >= RookValueMg;
}

bool is_KBPsK(const Position& pos, Color us) {
    return pos.non_pawn_material(us) == BishopValueMg && pos.count<PAWN>(us) >= 1;
}

bool is_KQKRPs(const Position& pos, Color us) {
    return !pos.count<PAWN>(us) && pos.non_pawn_material(us) == QueenValueMg
        && pos.count<ROOK>(~us) == 1 && pos.count<PAWN>(~us) >= 1;
}

template<Color Us>
EvalScore imbalance(const int pieceCount[][PIECE_TYPE_NB]) {

    constexpr Color Them = ~Us;

    EvalScore bonus = 0;

    for (int pt1 = NO_PIECE_TYPE; pt1 <= QUEEN; ++pt1)
    {
        if (!pieceCount[Us][pt1])
            continue;

        int v = QuadraticOurs[pt1][pt1] * pieceCount[Us][pt1];

        for (int pt2 = NO_PIECE_TYPE; pt2 < pt1; ++pt2)
            v += QuadraticOurs[pt1][pt2] * pieceCount[Us][pt2]
               + QuadraticTheirs[pt1][pt2] * pieceCount[Them][pt2];

        bonus += pieceCount[Us][pt1] * v;
    }

    return bonus;
}

}  // namespace

namespace Material {

Entry* probe(const Position& pos) {

    Key    key = pos.material_key();
    thread_local Table localMaterialTable;
    Entry* e   = localMaterialTable[key];

    if (e->key == key)
        return e;

    *e = Entry{};
    e->key = key;
    e->factor[WHITE] = e->factor[BLACK] = (uint8_t) SCALE_FACTOR_NORMAL;

    Value npm_w = pos.non_pawn_material(WHITE);
    Value npm_b = pos.non_pawn_material(BLACK);
    Value npm   = std::clamp(npm_w + npm_b, EndgameLimit, MidgameLimit);

    e->gamePhase = Phase(((npm - EndgameLimit) * PHASE_MIDGAME) / (MidgameLimit - EndgameLimit));

    if ((e->evaluationFunction = Endgames::probe<Value>(key)) != nullptr)
        return e;

    for (Color c : {WHITE, BLACK})
        if (is_KXK(pos, c))
        {
            e->evaluationFunction = &EvaluateKXK[c];
            return e;
        }

    const auto* sf = Endgames::probe<ScaleFactor>(key);

    if (sf)
    {
        e->scalingFunction[sf->strongSide] = sf;
        return e;
    }

    for (Color c : {WHITE, BLACK})
    {
        if (is_KBPsK(pos, c))
            e->scalingFunction[c] = &ScaleKBPsK[c];

        else if (is_KQKRPs(pos, c))
            e->scalingFunction[c] = &ScaleKQKRPs[c];
    }

    if (npm_w + npm_b == VALUE_ZERO && pos.pieces(PAWN))
    {
        if (!pos.count<PAWN>(BLACK))
        {
            assert(pos.count<PAWN>(WHITE) >= 2);
            e->scalingFunction[WHITE] = &ScaleKPsK[WHITE];
        }
        else if (!pos.count<PAWN>(WHITE))
        {
            assert(pos.count<PAWN>(BLACK) >= 2);
            e->scalingFunction[BLACK] = &ScaleKPsK[BLACK];
        }
        else if (pos.count<PAWN>(WHITE) == 1 && pos.count<PAWN>(BLACK) == 1)
        {
            e->scalingFunction[WHITE] = &ScaleKPKP[WHITE];
            e->scalingFunction[BLACK] = &ScaleKPKP[BLACK];
        }
    }

    if (!pos.count<PAWN>(WHITE) && npm_w - npm_b <= BishopValueMg)
        e->factor[WHITE] = uint8_t(npm_w < RookValueMg ? SCALE_FACTOR_DRAW
                                 : npm_b <= BishopValueMg ? 4 : 14);

    if (!pos.count<PAWN>(BLACK) && npm_b - npm_w <= BishopValueMg)
        e->factor[BLACK] = uint8_t(npm_b < RookValueMg ? SCALE_FACTOR_DRAW
                                 : npm_w <= BishopValueMg ? 4 : 14);

    const int pieceCount[COLOR_NB][PIECE_TYPE_NB] = {
      {pos.count<BISHOP>(WHITE) > 1, pos.count<PAWN>(WHITE), pos.count<KNIGHT>(WHITE),
       pos.count<BISHOP>(WHITE), pos.count<ROOK>(WHITE), pos.count<QUEEN>(WHITE)},
      {pos.count<BISHOP>(BLACK) > 1, pos.count<PAWN>(BLACK), pos.count<KNIGHT>(BLACK),
       pos.count<BISHOP>(BLACK), pos.count<ROOK>(BLACK), pos.count<QUEEN>(BLACK)}};

    e->score = (imbalance<WHITE>(pieceCount) - imbalance<BLACK>(pieceCount)) / 16;
    return e;
}

}  // namespace Material

}  // namespace Stockfish
