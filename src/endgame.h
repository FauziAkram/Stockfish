/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
*/

#ifndef ENDGAME_H_INCLUDED
#define ENDGAME_H_INCLUDED

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "position.h"
#include "types.h"

namespace Stockfish {

enum EndgameCode {

    EVALUATION_FUNCTIONS,
    KNNK,
    KNNKP,
    KXK,
    KBNK,
    KPK,
    KRKP,
    KRKB,
    KRKN,
    KQKP,
    KQKR,

    SCALING_FUNCTIONS,
    KBPsK,
    KQKRPs,
    KRPKR,
    KRPKB,
    KRPPKRP,
    KPsK,
    KBPKB,
    KBPPKB,
    KBPKN,
    KPKP
};

using EvalScore = int;
using ScaleFactor = int;

constexpr ScaleFactor SCALE_FACTOR_DRAW   = 0;
constexpr ScaleFactor SCALE_FACTOR_NORMAL = 64;
constexpr ScaleFactor SCALE_FACTOR_MAX    = 64;
constexpr ScaleFactor SCALE_FACTOR_NONE   = 255;

using Phase = int;
constexpr Phase PHASE_ENDGAME = 0;
constexpr Phase PHASE_MIDGAME = 128;

constexpr Value PawnValueEg   = PawnValue;
constexpr Value KnightValueMg = KnightValue;
constexpr Value BishopValueMg = BishopValue;
constexpr Value RookValueMg   = RookValue;
constexpr Value RookValueEg   = RookValue;
constexpr Value QueenValueMg  = QueenValue;
constexpr Value QueenValueEg  = QueenValue;
constexpr Value EndgameLimit  = 3915;
constexpr Value MidgameLimit  = 15258;
constexpr Value VALUE_KNOWN_WIN = 10000;

template<EndgameCode E>
using eg_type = typename std::conditional<(E < SCALING_FUNCTIONS), Value, ScaleFactor>::type;

template<typename T>
struct EndgameBase {

    explicit EndgameBase(Color c) : strongSide(c), weakSide(~c) {}
    virtual ~EndgameBase() = default;
    virtual T operator()(const Position&) const = 0;

    const Color strongSide, weakSide;
};

template<EndgameCode E, typename T = eg_type<E>>
struct Endgame : public EndgameBase<T> {

    explicit Endgame(Color c) : EndgameBase<T>(c) {}
    T operator()(const Position&) const override;
};

namespace Endgames {

template<typename T>
using Ptr = std::unique_ptr<EndgameBase<T>>;
template<typename T>
using Map = std::unordered_map<Key, Ptr<T>>;

extern std::pair<Map<Value>, Map<ScaleFactor>> maps;

void init();

template<typename T>
Map<T>& map() {
    return std::get<std::is_same<T, ScaleFactor>::value>(maps);
}

template<EndgameCode E, typename T = eg_type<E>>
void add(const std::string& code) {

    StateInfo st;
    map<T>()[Position().set(code, WHITE, &st).material_key()] = Ptr<T>(new Endgame<E>(WHITE));
    map<T>()[Position().set(code, BLACK, &st).material_key()] = Ptr<T>(new Endgame<E>(BLACK));
}

template<typename T>
const EndgameBase<T>* probe(Key key) {
    auto it = map<T>().find(key);
    return it != map<T>().end() ? it->second.get() : nullptr;
}

}  // namespace Endgames

namespace Bitbases {
bool probe(Square strongKing, Square strongPawn, Square weakKing, Color us);
}

constexpr EvalScore make_score(int mg, int eg) {
    return EvalScore((mg & 0xFFFF) + int(unsigned(eg & 0xFFFF) << 16));
}

constexpr int mg_value(EvalScore s) { return int16_t(s); }
constexpr int eg_value(EvalScore s) { return int16_t(unsigned(s + 0x8000) >> 16); }

}  // namespace Stockfish

#endif  // #ifndef ENDGAME_H_INCLUDED
