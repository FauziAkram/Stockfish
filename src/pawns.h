/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
*/

#ifndef PAWNS_H_INCLUDED
#define PAWNS_H_INCLUDED

#include <array>

#include "endgame.h"
#include "misc.h"
#include "position.h"
#include "types.h"

namespace Stockfish::Pawns {

struct Entry {

    EvalScore pawn_score(Color c) const { return scores[c]; }
    Bitboard  pawn_attacks(Color c) const { return pawnAttacks[c]; }
    Bitboard  passed_pawns(Color c) const { return passedPawns[c]; }
    Bitboard  pawn_attacks_span(Color c) const { return pawnAttacksSpan[c]; }
    int       passed_count() const { return popcount(passedPawns[WHITE] | passedPawns[BLACK]); }
    int       blocked_count() const { return blockedCount; }

    template<Color Us>
    EvalScore king_safety(const Position& pos) {
        int castling = int(pos.state()->castlingRights & (Us == WHITE ? WHITE_CASTLING : BLACK_CASTLING));
        return kingSquares[Us] == pos.square<KING>(Us) && castlingRights[Us] == castling
             ? kingSafety[Us]
             : (kingSafety[Us] = do_king_safety<Us>(pos));
    }

    template<Color Us>
    EvalScore do_king_safety(const Position& pos);

    template<Color Us>
    EvalScore evaluate_shelter(const Position& pos, Square ksq) const;

    Key      key = 0;
    EvalScore scores[COLOR_NB] = {0, 0};
    Bitboard passedPawns[COLOR_NB] = {0, 0};
    Bitboard pawnAttacks[COLOR_NB] = {0, 0};
    Bitboard pawnAttacksSpan[COLOR_NB] = {0, 0};
    Square   kingSquares[COLOR_NB] = {SQ_NONE, SQ_NONE};
    EvalScore kingSafety[COLOR_NB] = {0, 0};
    int      castlingRights[COLOR_NB] = {0, 0};
    int      blockedCount = 0;
};

class Table {
   public:
    Entry* operator[](Key key) { return &entries[key & (Size - 1)]; }

   private:
    static constexpr size_t         Size = 131072;
    std::array<Entry, Size> entries{};
};

Entry* probe(const Position& pos);

}  // namespace Stockfish::Pawns

#endif  // #ifndef PAWNS_H_INCLUDED
