/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
*/

#ifndef ENDGAME_H_INCLUDED
#define ENDGAME_H_INCLUDED

#include "types.h"

namespace Stockfish {

using EvalScore = int;
using ScaleFactor = int;

enum Phase : int { PHASE_ENDGAME = 0, PHASE_MIDGAME = 128 };

constexpr ScaleFactor SCALE_FACTOR_DRAW   = 0;
constexpr ScaleFactor SCALE_FACTOR_NORMAL = 64;
constexpr ScaleFactor SCALE_FACTOR_NONE   = 255;

constexpr Value BishopValueMg = BishopValue;
constexpr Value RookValueMg   = RookValue;
constexpr Value QueenValueMg  = QueenValue;
constexpr Value EndgameLimit  = 3915;
constexpr Value MidgameLimit  = 15258;

constexpr EvalScore make_score(int mg, int eg) { return (mg + eg) / 2; }

class Position;

template<typename T>
struct EndgameBase {
    explicit EndgameBase(Color c) : strongSide(c) {}
    virtual ~EndgameBase() = default;
    virtual T operator()(const Position&) const = 0;
    Color strongSide;
};

struct KXK {};
struct KBPsK {};
struct KQKRPs {};
struct KPsK {};
struct KPKP {};

template<typename Eg>
struct Endgame : EndgameBase<Value> {
    explicit Endgame(Color c) : EndgameBase<Value>(c) {}
    Value operator()(const Position&) const override { return VALUE_DRAW; }
};

template<>
struct Endgame<KBPsK> : EndgameBase<ScaleFactor> {
    explicit Endgame(Color c) : EndgameBase<ScaleFactor>(c) {}
    ScaleFactor operator()(const Position&) const override { return SCALE_FACTOR_NONE; }
};

template<>
struct Endgame<KQKRPs> : EndgameBase<ScaleFactor> {
    explicit Endgame(Color c) : EndgameBase<ScaleFactor>(c) {}
    ScaleFactor operator()(const Position&) const override { return SCALE_FACTOR_NONE; }
};

template<>
struct Endgame<KPsK> : EndgameBase<ScaleFactor> {
    explicit Endgame(Color c) : EndgameBase<ScaleFactor>(c) {}
    ScaleFactor operator()(const Position&) const override { return SCALE_FACTOR_NONE; }
};

template<>
struct Endgame<KPKP> : EndgameBase<ScaleFactor> {
    explicit Endgame(Color c) : EndgameBase<ScaleFactor>(c) {}
    ScaleFactor operator()(const Position&) const override { return SCALE_FACTOR_NONE; }
};

namespace Endgames {
template<typename T>
const EndgameBase<T>* probe(Key) {
    return nullptr;
}
}  // namespace Endgames

}  // namespace Stockfish

#endif
