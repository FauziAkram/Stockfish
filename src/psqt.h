/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2023 The Stockfish developers (see AUTHORS file)

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef PSQT_H_INCLUDED
#define PSQT_H_INCLUDED


#include "types.h"


namespace Stockfish::PSQT
{
  
// Parameters for Bishop tiling
extern const Score BishTiling[4];
// Parameters for Rook tiling
extern const Score RookTiling[4];
// Parameters for Queen tiling
extern const Score QueenTiling[4];

namespace hidden 
{
extern Score psq[PIECE_NB][SQUARE_NB];
  
  }

Piece idx(Piece p);

Score psq(Piece p, Square s);

// Fill psqt array from a set of internally linked parameters
void init();

} // namespace Stockfish::PSQT


#endif // PSQT_H_INCLUDED
