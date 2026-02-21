/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2026 The Stockfish developers (see AUTHORS file)
*/

#ifndef EVALUATE_H_INCLUDED
#define EVALUATE_H_INCLUDED

#include <string>

#include "types.h"

namespace Stockfish {

class Position;

namespace Eval {

std::string trace(Position& pos);
Value       evaluate(const Position& pos);

}  // namespace Eval

}  // namespace Stockfish

#endif  // #ifndef EVALUATE_H_INCLUDED
