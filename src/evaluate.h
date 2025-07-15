#ifndef EVALUATE_H_INCLUDED
#define EVALUATE_H_INCLUDED

#include <string>
#include "types.h"

namespace Stockfish {

class Position;

namespace Eval {

  Value evaluate(const Position& pos);
  std::string trace(Position& pos);

} // namespace Eval

} // namespace Stockfish

#endif // #ifndef EVALUATE_H_INCLUDED
