#ifndef KPK_BITBASE_H_INCLUDED
#define KPK_BITBASE_H_INCLUDED

#include "types.h"

namespace Stockfish {

namespace KpkBitbase {

// Initializes the KPK bitbase. This should be called once at startup.
void init();

// Probes the KPK bitbase for a win for White.
// IMPORTANT: Assumes the pawn 'wpsq' belongs to White.
// Returns true if the position is a win for White, false otherwise (draw or loss).
bool probe(Square wksq, Square wpsq, Square bksq, Color stm);

} // namespace KpkBitbase

} // namespace Stockfish

#endif // KPK_BITBASE_H_INCLUDED
