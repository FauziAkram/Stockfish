/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
*/

#ifndef MATERIAL_H_INCLUDED
#define MATERIAL_H_INCLUDED

#include <array>

#include "endgame.h"
#include "position.h"
#include "types.h"

namespace Stockfish::Material {

struct Entry {

    EvalScore imbalance() const { return score; }
    Phase     game_phase() const { return (Phase) gamePhase; }
    bool      specialized_eval_exists() const { return evaluationFunction != nullptr; }
    Value     evaluate(const Position& pos) const { return (*evaluationFunction)(pos); }

    ScaleFactor scale_factor(const Position& pos, Color c) const {
        ScaleFactor sf = scalingFunction[c] ? (*scalingFunction[c])(pos) : SCALE_FACTOR_NONE;
        return sf != SCALE_FACTOR_NONE ? sf : ScaleFactor(factor[c]);
    }

    Key                           key = 0;
    const EndgameBase<Value>*     evaluationFunction = nullptr;
    const EndgameBase<ScaleFactor>* scalingFunction[COLOR_NB] = {nullptr, nullptr};
    EvalScore                     score = 0;
    int16_t                       gamePhase = 0;
    uint8_t                       factor[COLOR_NB] = {SCALE_FACTOR_NORMAL, SCALE_FACTOR_NORMAL};
};

class Table {
   public:
    Entry* operator[](Key key) { return &entries[key & (Size - 1)]; }

   private:
    static constexpr size_t         Size = 8192;
    std::array<Entry, Size> entries{};
};

Entry* probe(const Position& pos);

}  // namespace Stockfish::Material

#endif  // #ifndef MATERIAL_H_INCLUDED
