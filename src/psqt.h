#ifndef PSQT_H_INCLUDED
#define PSQT_H_INCLUDED

#include "types.h"

namespace Stockfish::PSQT {

extern Score psq[PIECE_NB][SQUARE_NB];
void init();

} // namespace Stockfish::PSQT

#endif // PSQT_H_INCLUDED
