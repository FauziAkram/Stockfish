#ifndef BITBASES_H_INCLUDED
#define BITBASES_H_INCLUDED

#include "types.h"

namespace Stockfish {

namespace Bitbases {

void init();
bool probe(Square wksq, Square wpsq, Square bksq, Color us);

} // namespace Bitbases

} // namespace Stockfish

#endif // #ifndef BITBASES_H_INCLUDED
