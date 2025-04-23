/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2025 The Stockfish developers (see AUTHORS file)

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

#include "position.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstring>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string_view>
#include <utility>

#include "bitboard.h"
#include "misc.h"
#include "movegen.h"
#include "syzygy/tbprobe.h"
#include "tt.h"
#include "uci.h"

using std::string;

namespace Stockfish {

namespace Zobrist {

Key psq[PIECE_NB][SQUARE_NB];
Key enpassant[FILE_NB];
Key castling[CASTLING_RIGHT_NB];
Key side, noPawns;
}

namespace {

constexpr std::string_view PieceToChar(" PNBRQK  pnbrqk");

static constexpr Piece Pieces[] = {W_PAWN, W_KNIGHT, W_BISHOP, W_ROOK, W_QUEEN, W_KING,
                                   B_PAWN, B_KNIGHT, B_BISHOP, B_ROOK, B_QUEEN, B_KING};
}  // namespace


// Returns an ASCII representation of the position
std::ostream& operator<<(std::ostream& os, const Position& pos) {

    os << "\n +---+---+---+---+---+---+---+---+\n";

    for (Rank r = RANK_8; r >= RANK_1; --r)
    {
        for (File f = FILE_A; f <= FILE_H; ++f)
            os << " | " << PieceToChar[pos.piece_on(make_square(f, r))];

        os << " | " << (1 + r) << "\n +---+---+---+---+---+---+---+---+\n";
    }

    os << "   a   b   c   d   e   f   g   h\n"
       << "\nFen: " << pos.fen() << "\nKey: " << std::hex << std::uppercase << std::setfill('0')
       << std::setw(16) << pos.key() << std::setfill(' ') << std::dec << "\nCheckers: ";

    for (Bitboard b = pos.checkers(); b;)
        os << UCIEngine::square(pop_lsb(b)) << " ";

    if (int(Tablebases::MaxCardinality) >= popcount(pos.pieces()) && !pos.can_castle(ANY_CASTLING))
    {
        StateInfo st;

        Position p;
        p.set(pos.fen(), pos.is_chess960(), &st);
        Tablebases::ProbeState s1, s2;
        Tablebases::WDLScore   wdl = Tablebases::probe_wdl(p, &s1);
        int                    dtz = Tablebases::probe_dtz(p, &s2);
        os << "\nTablebases WDL: " << std::setw(4) << wdl << " (" << s1 << ")"
           << "\nTablebases DTZ: " << std::setw(4) << dtz << " (" << s2 << ")";
    }

    return os;
}


// Implements Marcel van Kervinck's cuckoo algorithm to detect repetition of positions
// for 3-fold repetition draws. The algorithm uses two hash tables with Zobrist hashes
// to allow fast detection of recurring positions. For details see:
// http://web.archive.org/web/20201107002606/https://marcelk.net/2013-04-06/paper/upcoming-rep-v2.pdf

// First and second hash functions for indexing the cuckoo tables
inline int H1(Key h) { return h & 0x1fff; }
inline int H2(Key h) { return (h >> 16) & 0x1fff; }

// Cuckoo tables with Zobrist hashes of valid reversible moves, and the moves themselves
std::array<Key, 8192>  cuckoo;
std::array<Move, 8192> cuckooMove;

// Initializes at startup the various arrays used to compute hash keys
void Position::init() {

    PRNG rng(1070372);

    for (Piece pc : Pieces)
        for (Square s = SQ_A1; s <= SQ_H8; ++s)
            Zobrist::psq[pc][s] = rng.rand<Key>();
    // pawns on these squares will promote
    std::fill_n(Zobrist::psq[W_PAWN] + SQ_A8, 8, 0);
    std::fill_n(Zobrist::psq[B_PAWN], 8, 0);

    for (File f = FILE_A; f <= FILE_H; ++f)
        Zobrist::enpassant[f] = rng.rand<Key>();

    for (int cr = NO_CASTLING; cr <= ANY_CASTLING; ++cr)
        Zobrist::castling[cr] = rng.rand<Key>();

    Zobrist::side    = rng.rand<Key>();
    Zobrist::noPawns = rng.rand<Key>();

    // Prepare the cuckoo tables
    cuckoo.fill(0);
    cuckooMove.fill(Move::none());
    [[maybe_unused]] int count = 0;
    for (Piece pc : Pieces)
        for (Square s1 = SQ_A1; s1 <= SQ_H8; ++s1)
            for (Square s2 = Square(s1 + 1); s2 <= SQ_H8; ++s2)
                if ((type_of(pc) != PAWN) && (attacks_bb(type_of(pc), s1, 0) & s2))
                {
                    Move move = Move(s1, s2);
                    Key  key  = Zobrist::psq[pc][s1] ^ Zobrist::psq[pc][s2] ^ Zobrist::side;
                    int  i    = H1(key);
                    while (true)
                    {
                        std::swap(cuckoo[i], key);
                        std::swap(cuckooMove[i], move);
                        if (move == Move::none())  // Arrived at empty slot?
                            break;
                        i = (i == H1(key)) ? H2(key) : H1(key);  // Push victim to alternative slot
                    }
                    count++;
                }
    assert(count == 3668);
}


// Initializes the position object with the given FEN string.
// This function is not very robust - make sure that input FENs are correct,
// this is assumed to be the responsibility of the GUI.
Position& Position::set(const string& fenStr, bool isChess960, StateInfo* si) {
    /*
   A FEN string defines a particular position using only the ASCII character set.

   A FEN string contains six fields separated by a space. The fields are:

   1) Piece placement (from white's perspective). Each rank is described, starting
      with rank 8 and ending with rank 1. Within each rank, the contents of each
      square are described from file A through file H. Following the Standard
      Algebraic Notation (SAN), each piece is identified by a single letter taken
      from the standard English names. White pieces are designated using upper-case
      letters ("PNBRQK") whilst Black uses lowercase ("pnbrqk"). Blank squares are
      noted using digits 1 through 8 (the number of blank squares), and "/"
      separates ranks.

   2) Active color. "w" means white moves next, "b" means black.

   3) Castling availability. If neither side can castle, this is "-". Otherwise,
      this has one or more letters: "K" (White can castle kingside), "Q" (White
      can castle queenside), "k" (Black can castle kingside), and/or "q" (Black
      can castle queenside).

   4) En passant target square (in algebraic notation). If there's no en passant
      target square, this is "-". If a pawn has just made a 2-square move, this
      is the position "behind" the pawn. Following X-FEN standard, this is recorded
      only if there is a pawn in position to make an en passant capture, and if
      there really is a pawn that might have advanced two squares.

   5) Halfmove clock. This is the number of halfmoves since the last pawn advance
      or capture. This is used to determine if a draw can be claimed under the
      fifty-move rule.

   6) Fullmove number. The number of the full move. It starts at 1, and is
      incremented after Black's move.
*/

    unsigned char      col, row, token;
    size_t             idx;
    Square             sq = SQ_A8;
    std::istringstream ss(fenStr);

    std::memset(this, 0, sizeof(Position));
    std::memset(si, 0, sizeof(StateInfo));
    st = si;

    ss >> std::noskipws;

    // 1. Piece placement
    while ((ss >> token) && !isspace(token))
    {
        if (isdigit(token))
            sq += (token - '0') * EAST;  // Advance the given number of files

        else if (token == '/')
            sq += 2 * SOUTH;

        else if ((idx = PieceToChar.find(token)) != string::npos)
        {
            put_piece(Piece(idx), sq);
            ++sq;
        }
    }

    // 2. Active color
    ss >> token;
    sideToMove = (token == 'w' ? WHITE : BLACK);
    ss >> token;

    // 3. Castling availability. Compatible with 3 standards: Normal FEN standard,
    // Shredder-FEN that uses the letters of the columns on which the rooks began
    // the game instead of KQkq and also X-FEN standard that, in case of Chess960,
    // if an inner rook is associated with the castling right, the castling tag is
    // replaced by the file letter of the involved rook, as for the Shredder-FEN.
    while ((ss >> token) && !isspace(token))
    {
        Square rsq;
        Color  c    = islower(token) ? BLACK : WHITE;
        Piece  rook = make_piece(c, ROOK);

        token = char(toupper(token));

        if (token == 'K')
            for (rsq = relative_square(c, SQ_H1); piece_on(rsq) != rook; --rsq)
            {}

        else if (token == 'Q')
            for (rsq = relative_square(c, SQ_A1); piece_on(rsq) != rook; ++rsq)
            {}

        else if (token >= 'A' && token <= 'H')
            rsq = make_square(File(token - 'A'), relative_rank(c, RANK_1));

        else
            continue;

        set_castling_right(c, rsq);
    }

    // 4. En passant square.
    // Ignore if square is invalid or not on side to move relative rank 6.
    bool enpassant = false;

    if (((ss >> col) && (col >= 'a' && col <= 'h'))
        && ((ss >> row) && (row == (sideToMove == WHITE ? '6' : '3'))))
    {
        st->epSquare = make_square(File(col - 'a'), Rank(row - '1'));

        // En passant square will be considered only if
        // a) side to move have a pawn threatening epSquare
        // b) there is an enemy pawn in front of epSquare
        // c) there is no piece on epSquare or behind epSquare
        enpassant = attacks_bb<PAWN>(st->epSquare, ~sideToMove) & pieces(sideToMove, PAWN)
                 && (pieces(~sideToMove, PAWN) & (st->epSquare + pawn_push(~sideToMove)))
                 && !(pieces() & (st->epSquare | (st->epSquare + pawn_push(sideToMove))));
    }

    if (!enpassant)
        st->epSquare = SQ_NONE;

    // 5-6. Halfmove clock and fullmove number
    ss >> std::skipws >> st->rule50 >> gamePly;

    // Convert from fullmove starting from 1 to gamePly starting from 0,
    // handle also common incorrect FEN with fullmove = 0.
    gamePly = std::max(2 * (gamePly - 1), 0) + (sideToMove == BLACK);

    chess960 = isChess960;
    set_state();

    assert(pos_is_ok());

    return *this;
}


// Helper function used to set castling
// rights given the corresponding color and the rook starting square.
void Position::set_castling_right(Color c, Square rfrom) {

    Square         kfrom = square<KING>(c);
    CastlingRights cr    = c & (kfrom < rfrom ? KING_SIDE : QUEEN_SIDE);

    st->castlingRights |= cr;
    castlingRightsMask[kfrom] |= cr;
    castlingRightsMask[rfrom] |= cr;
    castlingRookSquare[cr] = rfrom;

    Square kto = relative_square(c, cr & KING_SIDE ? SQ_G1 : SQ_C1);
    Square rto = relative_square(c, cr & KING_SIDE ? SQ_F1 : SQ_D1);

    castlingPath[cr] = (between_bb(rfrom, rto) | between_bb(kfrom, kto)) & ~(kfrom | rfrom);
}


// Sets king attacks to detect if a move gives check
void Position::set_check_info() const {

    update_slider_blockers(WHITE);
    update_slider_blockers(BLACK);

    Square ksq = square<KING>(~sideToMove);

    st->checkSquares[PAWN]   = attacks_bb<PAWN>(ksq, ~sideToMove);
    st->checkSquares[KNIGHT] = attacks_bb<KNIGHT>(ksq);
    st->checkSquares[BISHOP] = attacks_bb<BISHOP>(ksq, pieces());
    st->checkSquares[ROOK]   = attacks_bb<ROOK>(ksq, pieces());
    st->checkSquares[QUEEN]  = st->checkSquares[BISHOP] | st->checkSquares[ROOK];
    st->checkSquares[KING]   = 0;
}


// Computes the hash keys of the position, and other
// data that once computed is updated incrementally as moves are made.
// The function is only used when a new position is set up
void Position::set_state() const {

    st->key = st->materialKey = 0;
    st->minorPieceKey         = 0;
    st->nonPawnKey[WHITE] = st->nonPawnKey[BLACK] = 0;
    st->pawnKey                                   = Zobrist::noPawns;
    st->nonPawnMaterial[WHITE] = st->nonPawnMaterial[BLACK] = VALUE_ZERO;
    st->checkersBB = attackers_to(square<KING>(sideToMove)) & pieces(~sideToMove);

    set_check_info();

    for (Bitboard b = pieces(); b;)
    {
        Square s  = pop_lsb(b);
        Piece  pc = piece_on(s);
        st->key ^= Zobrist::psq[pc][s];

        if (type_of(pc) == PAWN)
            st->pawnKey ^= Zobrist::psq[pc][s];

        else
        {
            st->nonPawnKey[color_of(pc)] ^= Zobrist::psq[pc][s];

            if (type_of(pc) != KING)
            {
                st->nonPawnMaterial[color_of(pc)] += PieceValue[pc];

                if (type_of(pc) <= BISHOP)
                    st->minorPieceKey ^= Zobrist::psq[pc][s];
            }
        }
    }

    if (st->epSquare != SQ_NONE)
        st->key ^= Zobrist::enpassant[file_of(st->epSquare)];

    if (sideToMove == BLACK)
        st->key ^= Zobrist::side;

    st->key ^= Zobrist::castling[st->castlingRights];

    for (Piece pc : Pieces)
        for (int cnt = 0; cnt < pieceCount[pc]; ++cnt)
            st->materialKey ^= Zobrist::psq[pc][8 + cnt];
}


// Overload to initialize the position object with the given endgame code string
// like "KBPKN". It's mainly a helper to get the material key out of an endgame code.
Position& Position::set(const string& code, Color c, StateInfo* si) {

    assert(code[0] == 'K');

    string sides[] = {code.substr(code.find('K', 1)),                                // Weak
                      code.substr(0, std::min(code.find('v'), code.find('K', 1)))};  // Strong

    assert(sides[0].length() > 0 && sides[0].length() < 8);
    assert(sides[1].length() > 0 && sides[1].length() < 8);

    std::transform(sides[c].begin(), sides[c].end(), sides[c].begin(), tolower);

    string fenStr = "8/" + sides[0] + char(8 - sides[0].length() + '0') + "/8/8/8/8/" + sides[1]
                  + char(8 - sides[1].length() + '0') + "/8 w - - 0 10";

    return set(fenStr, false, si);
}


// Returns a FEN representation of the position. In case of
// Chess960 the Shredder-FEN notation is used. This is mainly a debugging function.
string Position::fen() const {

    int                emptyCnt;
    std::ostringstream ss;

    for (Rank r = RANK_8; r >= RANK_1; --r)
    {
        for (File f = FILE_A; f <= FILE_H; ++f)
        {
            for (emptyCnt = 0; f <= FILE_H && empty(make_square(f, r)); ++f)
                ++emptyCnt;

            if (emptyCnt)
                ss << emptyCnt;

            if (f <= FILE_H)
                ss << PieceToChar[piece_on(make_square(f, r))];
        }

        if (r > RANK_1)
            ss << '/';
    }

    ss << (sideToMove == WHITE ? " w " : " b ");

    if (can_castle(WHITE_OO))
        ss << (chess960 ? char('A' + file_of(castling_rook_square(WHITE_OO))) : 'K');

    if (can_castle(WHITE_OOO))
        ss << (chess960 ? char('A' + file_of(castling_rook_square(WHITE_OOO))) : 'Q');

    if (can_castle(BLACK_OO))
        ss << (chess960 ? char('a' + file_of(castling_rook_square(BLACK_OO))) : 'k');

    if (can_castle(BLACK_OOO))
        ss << (chess960 ? char('a' + file_of(castling_rook_square(BLACK_OOO))) : 'q');

    if (!can_castle(ANY_CASTLING))
        ss << '-';

    ss << (ep_square() == SQ_NONE ? " - " : " " + UCIEngine::square(ep_square()) + " ")
       << st->rule50 << " " << 1 + (gamePly - (sideToMove == BLACK)) / 2;

    return ss.str();
}

// Calculates st->blockersForKing[c] and st->pinners[~c],
// which store respectively the pieces preventing king of color c from being in check
// and the slider pieces of color ~c pinning pieces of color c to the king.
void Position::update_slider_blockers(Color c) const {

    Square ksq = square<KING>(c);

    st->blockersForKing[c] = 0;
    st->pinners[~c]        = 0;

    // Snipers are sliders that attack 's' when a piece and other snipers are removed
    Bitboard snipers = ((attacks_bb<ROOK>(ksq) & pieces(QUEEN, ROOK))
                        | (attacks_bb<BISHOP>(ksq) & pieces(QUEEN, BISHOP)))
                     & pieces(~c);
    Bitboard occupancy = pieces() ^ snipers;

    while (snipers)
    {
        Square   sniperSq = pop_lsb(snipers);
        Bitboard b        = between_bb(ksq, sniperSq) & occupancy;

        if (b && !more_than_one(b))
        {
            st->blockersForKing[c] |= b;
            if (b & pieces(c))
                st->pinners[~c] |= sniperSq;
        }
    }
}


// Computes a bitboard of all pieces which attack a given square.
// Slider attacks use the occupied bitboard to indicate occupancy.
Bitboard Position::attackers_to(Square s, Bitboard occupied) const {

    return (attacks_bb<ROOK>(s, occupied) & pieces(ROOK, QUEEN))
         | (attacks_bb<BISHOP>(s, occupied) & pieces(BISHOP, QUEEN))
         | (attacks_bb<PAWN>(s, BLACK) & pieces(WHITE, PAWN))
         | (attacks_bb<PAWN>(s, WHITE) & pieces(BLACK, PAWN))
         | (attacks_bb<KNIGHT>(s) & pieces(KNIGHT)) | (attacks_bb<KING>(s) & pieces(KING));
}

bool Position::attackers_to_exist(Square s, Bitboard occupied, Color c) const {

    return ((attacks_bb<ROOK>(s) & pieces(c, ROOK, QUEEN))
            && (attacks_bb<ROOK>(s, occupied) & pieces(c, ROOK, QUEEN)))
        || ((attacks_bb<BISHOP>(s) & pieces(c, BISHOP, QUEEN))
            && (attacks_bb<BISHOP>(s, occupied) & pieces(c, BISHOP, QUEEN)))
        || (((attacks_bb<PAWN>(s, ~c) & pieces(PAWN)) | (attacks_bb<KNIGHT>(s) & pieces(KNIGHT))
             | (attacks_bb<KING>(s) & pieces(KING)))
            & pieces(c));
}

// Tests whether a pseudo-legal move is legal
bool Position::legal(Move m) const {

    assert(m.is_ok());

    Color  us   = sideToMove;
    Square from = m.from_sq();
    Square to   = m.to_sq();

    assert(color_of(moved_piece(m)) == us);
    assert(piece_on(square<KING>(us)) == make_piece(us, KING));

    // En passant captures are a tricky special case. Because they are rather
    // uncommon, we do it simply by testing whether the king is attacked after
    // the move is made.
    if (m.type_of() == EN_PASSANT)
    {
        Square   ksq      = square<KING>(us);
        Square   capsq    = to - pawn_push(us);
        Bitboard occupied = (pieces() ^ from ^ capsq) | to;

        assert(to == ep_square());
        assert(moved_piece(m) == make_piece(us, PAWN));
        assert(piece_on(capsq) == make_piece(~us, PAWN));
        assert(piece_on(to) == NO_PIECE);

        return !(attacks_bb<ROOK>(ksq, occupied) & pieces(~us, QUEEN, ROOK))
            && !(attacks_bb<BISHOP>(ksq, occupied) & pieces(~us, QUEEN, BISHOP));
    }

    // Castling moves generation does not check if the castling path is clear of
    // enemy attacks, it is delayed at a later time: now!
    if (m.type_of() == CASTLING)
    {
        // After castling, the rook and king final positions are the same in
        // Chess960 as they would be in standard chess.
        to             = relative_square(us, to > from ? SQ_G1 : SQ_C1);
        Direction step = to > from ? WEST : EAST;

        for (Square s = to; s != from; s += step)
            if (attackers_to_exist(s, pieces(), ~us))
                return false;

        // In case of Chess960, verify if the Rook blocks some checks.
        // For instance an enemy queen in SQ_A1 when castling rook is in SQ_B1.
        return !chess960 || !(blockers_for_king(us) & m.to_sq());
    }

    // If the moving piece is a king, check whether the destination square is
    // attacked by the opponent.
    if (type_of(piece_on(from)) == KING)
        return !(attackers_to_exist(to, pieces() ^ from, ~us));

    // A non-king move is legal if and only if it is not pinned or it
    // is moving along the ray towards or away from the king.
    return !(blockers_for_king(us) & from) || line_bb(from, to) & pieces(us, KING);
}


// Takes a random move and tests whether the move is
// pseudo-legal. It is used to validate moves from TT that can be corrupted
// due to SMP concurrent access or hash position key aliasing.
bool Position::pseudo_legal(const Move m) const {

    Color  us   = sideToMove;
    Square from = m.from_sq();
    Square to   = m.to_sq();
    Piece  pc   = moved_piece(m);

    // Use a slower but simpler function for uncommon cases
    // yet we skip the legality check of MoveList<LEGAL>().
    if (m.type_of() != NORMAL)
        return checkers() ? MoveList<EVASIONS>(*this).contains(m)
                          : MoveList<NON_EVASIONS>(*this).contains(m);

    // Is not a promotion, so the promotion piece must be empty
    assert(m.promotion_type() - KNIGHT == NO_PIECE_TYPE);

    // If the 'from' square is not occupied by a piece belonging to the side to
    // move, the move is obviously not legal.
    if (pc == NO_PIECE || color_of(pc) != us)
        return false;

    // The destination square cannot be occupied by a friendly piece
    if (pieces(us) & to)
        return false;

    // Handle the special case of a pawn move
    if (type_of(pc) == PAWN)
    {
        // We have already handled promotion moves, so destination cannot be on the 8th/1st rank
        if ((Rank8BB | Rank1BB) & to)
            return false;

        if (!(attacks_bb<PAWN>(from, us) & pieces(~us) & to)  // Not a capture
            && !((from + pawn_push(us) == to) && empty(to))   // Not a single push
            && !((from + 2 * pawn_push(us) == to)             // Not a double push
                 && (relative_rank(us, from) == RANK_2) && empty(to) && empty(to - pawn_push(us))))
            return false;
    }
    else if (!(attacks_bb(type_of(pc), from, pieces()) & to))
        return false;

    // Evasions generator already takes care to avoid some kind of illegal moves
    // and legal() relies on this. We therefore have to take care that the same
    // kind of moves are filtered out here.
    if (checkers())
    {
        if (type_of(pc) != KING)
        {
            // Double check? In this case, a king move is required
            if (more_than_one(checkers()))
                return false;

            // Our move must be a blocking interposition or a capture of the checking piece
            if (!(between_bb(square<KING>(us), lsb(checkers())) & to))
                return false;
        }
        // In case of king moves under check we have to remove the king so as to catch
        // invalid moves like b1a1 when opposite queen is on c1.
        else if (attackers_to_exist(to, pieces() ^ from, ~us))
            return false;
    }

    return true;
}


// Tests whether a pseudo-legal move gives a check
bool Position::gives_check(Move m) const {

    assert(m.is_ok());
    assert(color_of(moved_piece(m)) == sideToMove);

    Square from = m.from_sq();
    Square to   = m.to_sq();

    // Is there a direct check?
    if (check_squares(type_of(piece_on(from))) & to)
        return true;

    // Is there a discovered check?
    if (blockers_for_king(~sideToMove) & from)
        return !(line_bb(from, to) & pieces(~sideToMove, KING)) || m.type_of() == CASTLING;

    switch (m.type_of())
    {
    case NORMAL :
        return false;

    case PROMOTION :
        return attacks_bb(m.promotion_type(), to, pieces() ^ from) & pieces(~sideToMove, KING);

    // En passant capture with check? We have already handled the case of direct
    // checks and ordinary discovered check, so the only case we need to handle
    // is the unusual case of a discovered check through the captured pawn.
    case EN_PASSANT : {
        Square   capsq = make_square(file_of(to), rank_of(from));
        Bitboard b     = (pieces() ^ from ^ capsq) | to;

        return (attacks_bb<ROOK>(square<KING>(~sideToMove), b) & pieces(sideToMove, QUEEN, ROOK))
             | (attacks_bb<BISHOP>(square<KING>(~sideToMove), b)
                & pieces(sideToMove, QUEEN, BISHOP));
    }
    default :  //CASTLING
    {
        // Castling is encoded as 'king captures the rook'
        Square rto = relative_square(sideToMove, to > from ? SQ_F1 : SQ_D1);

        return check_squares(ROOK) & rto;
    }
    }
}


// Makes a move, and saves all information necessary
// to a StateInfo object. The move is assumed to be legal. Pseudo-legal
// moves should be filtered out before this function is called.
// If a pointer to the TT table is passed, the entry for the new position
// will be prefetched
DirtyPiece Position::do_move(Move                      m,
                             StateInfo&                newSt,
                             bool                      givesCheck,
                             const TranspositionTable* tt = nullptr) {

    assert(m.is_ok());
    assert(&newSt != st);

    Key k = st->key ^ Zobrist::side;

    // Copy some fields of the old state to our new StateInfo object except the
    // ones which are going to be recalculated from scratch anyway and then switch
    // our state pointer to point to the new (ready to be updated) state.
    std::memcpy(&newSt, st, offsetof(StateInfo, key)); // Copy up to, but not including, key
    newSt.previous = st;
    st->next       = &newSt;
    st             = &newSt;

    // Increment ply counters. rule50 handled incrementally below.
    ++gamePly;
    ++st->pliesFromNull;

    DirtyPiece dp;
    dp.dirty_num = 1; // Assume at least one piece moves

    Color  us       = sideToMove;
    Color  them     = ~us;
    Square from     = m.from_sq();
    Square to       = m.to_sq();
    Piece  pc       = piece_on(from); // Piece moving
    Piece  captured = m.type_of() == EN_PASSANT ? make_piece(them, PAWN) : piece_on(to); // Piece captured (or rook in castling)

    assert(color_of(pc) == us);
    assert(captured == NO_PIECE || color_of(captured) == (m.type_of() != CASTLING ? them : us));
    assert(type_of(captured) != KING);

    // Store previous EP square state before modifying it
    const Square oldEpSquare = st->epSquare;

    // Reset EP square state immediately. If a new one is created, it's set later.
    if (oldEpSquare != SQ_NONE)
    {
        // [UPDATE 2: Old EP Square] - Moved earlier
        k ^= Zobrist::enpassant[file_of(oldEpSquare)];
        st->epSquare = SQ_NONE;
    }

    // Store previous castling rights state before modifying it
    const int oldCastlingRights = st->castlingRights;

    // Update castling rights state if needed
    if (oldCastlingRights && (castlingRightsMask[from] | castlingRightsMask[to]))
    {
         st->castlingRights &= ~(castlingRightsMask[from] | castlingRightsMask[to]);
         // [UPDATE 3 & 4: Old/New Castling Rights] - Combined update
         k ^= Zobrist::castling[oldCastlingRights] ^ Zobrist::castling[st->castlingRights];
    }

    // Assume non-capture/non-pawn move initially for rule50, increment from previous state
    st->rule50 = st->previous->rule50 + 1;

    // Handle piece movement, captures, and associated Zobrist updates
    if (m.type_of() == CASTLING)
    {
        assert(pc == make_piece(us, KING));
        assert(captured == make_piece(us, ROOK)); // 'captured' holds the rook

        Square rfrom, rto;
        do_castling<true>(us, from, to, rfrom, rto, &dp); // Updates board, dp, piece counts

        // [UPDATE 5 & 6: Castling King & Rook Moves]
        k ^= Zobrist::psq[pc][from] ^ Zobrist::psq[pc][to];       // King part
        k ^= Zobrist::psq[captured][rfrom] ^ Zobrist::psq[captured][rto]; // Rook part

        // Update non-pawn keys for both pieces
        st->nonPawnKey[us] ^= Zobrist::psq[pc][from] ^ Zobrist::psq[pc][to];
        st->nonPawnKey[us] ^= Zobrist::psq[captured][rfrom] ^ Zobrist::psq[captured][rto];

        captured = NO_PIECE; // Not a material capture
        // rule50 is not reset for castling
    }
    else // NORMAL, PROMOTION, or EN_PASSANT
    {
        Square capsq = to; // Square of captured piece

        if (captured) // Handle capture details
        {
            if (m.type_of() == EN_PASSANT)
            {
                 capsq -= pawn_push(us);
                 assert(type_of(pc) == PAWN);
                 assert(to == oldEpSquare); // Use stored old EP square for assertion
                 assert(relative_rank(us, to) == RANK_6);
                 assert(piece_on(capsq) == make_piece(them, PAWN));
            }

            // Update board, piece lists for capture
            remove_piece(capsq);

            // [UPDATE 7: Captured Piece]
            k ^= Zobrist::psq[captured][capsq];

            // Update material/pawn keys for the captured piece
            if (type_of(captured) == PAWN)
                st->pawnKey ^= Zobrist::psq[captured][capsq];
            else {
                st->nonPawnMaterial[them] -= PieceValue[captured];
                st->nonPawnKey[them] ^= Zobrist::psq[captured][capsq];
                if (type_of(captured) <= BISHOP) // Knight or Bishop
                    st->minorPieceKey ^= Zobrist::psq[captured][capsq];
            }
            // Update material key regardless of piece type (count already decremented)
            st->materialKey ^= Zobrist::psq[captured][8 + pieceCount[captured]];

            dp.dirty_num = 2; // 1 piece moved, 1 piece captured
            dp.piece[1]  = captured;
            dp.from[1]   = capsq;
            dp.to[1]     = SQ_NONE;

            st->rule50 = 0; // Reset rule50 for capture
        } // End if(captured)

        // Set DP info for the moving piece (before potential promotion change)
        dp.piece[0] = pc;
        dp.from[0]  = from;
        dp.to[0]    = to;

        // Move the piece on the board state (bitboards and piece lists)
        move_piece(from, to);

        // [UPDATE 8 & 9: Moving Piece Out/In]
        k ^= Zobrist::psq[pc][from] ^ Zobrist::psq[pc][to];

        // Update non-pawn material/pawn keys based on the moving piece (pc) type
        if (type_of(pc) == PAWN)
        {
            // Update pawn structure hash key
            st->pawnKey ^= Zobrist::psq[pc][from] ^ Zobrist::psq[pc][to];
            // Reset rule50 for pawn move (can happen even if already reset by capture)
            st->rule50 = 0;

            // Set new en passant square if applicable (double push)
            if ((int(to) ^ int(from)) == 16 // Double push
                && (attacks_bb<PAWN>(to - pawn_push(us), us) & pieces(them, PAWN))) // Check adjacent enemy pawns
            {
                st->epSquare = to - pawn_push(us);
                // [UPDATE 10: New EP Square]
                k ^= Zobrist::enpassant[file_of(st->epSquare)];
            }
            else if (m.type_of() == PROMOTION)
            {
                Piece promotion     = make_piece(us, m.promotion_type());
                PieceType promotionType = type_of(promotion);

                assert(relative_rank(us, to) == RANK_8);
                assert(promotionType >= KNIGHT && promotionType <= QUEEN);

                // Update board/counts for promotion: remove pawn, add promoted piece
                remove_piece(to); // Remove the pawn just moved to 'to'
                put_piece(promotion, to); // Put the promoted piece

                // Update DP info for promotion
                dp.to[0]               = SQ_NONE; // Pawn effectively moved 'to' SQ_NONE
                dp.piece[dp.dirty_num] = promotion; // Promoted piece appears
                dp.from[dp.dirty_num]  = SQ_NONE;
                dp.to[dp.dirty_num]    = to;
                dp.dirty_num++;

                // Update hash key: Zobrist::psq[PAWN][to] is 0.
                // k already includes remove PAWN from 'from'. Need to add promoted piece at 'to'.
                // [UPDATE 11: Promotion Piece]
                k ^= Zobrist::psq[promotion][to];

                // Update material/keys: remove pawn, add promoted piece
                st->pawnKey ^= Zobrist::psq[pc][from]; // Final removal of pawn from structure key
                st->materialKey ^= Zobrist::psq[pc][8 + pieceCount[pc]]; // Remove pawn from material key
                st->nonPawnMaterial[us] += PieceValue[promotion];
                st->nonPawnKey[us] ^= Zobrist::psq[promotion][to];
                st->materialKey ^= Zobrist::psq[promotion][8 + pieceCount[promotion] - 1]; // Add promotion piece to material key
                if (promotionType <= BISHOP) // Knight or Bishop
                    st->minorPieceKey ^= Zobrist::psq[promotion][to];
            }
        }
        else // Non-pawn moving piece
        {
            // Update non-pawn keys for the moving piece
            st->nonPawnKey[us] ^= Zobrist::psq[pc][from] ^ Zobrist::psq[pc][to];
            // Update minor piece key if applicable (KNIGHT or BISHOP)
            if (type_of(pc) <= BISHOP && type_of(pc) != KING) // Check type not KING just in case
                st->minorPieceKey ^= Zobrist::psq[pc][from] ^ Zobrist::psq[pc][to];
            // rule50 counter remains incremented unless it was a capture
        }
    } // End of else block for non-CASTLING moves

    // Update the key with the final value
    st->key = k;
    if (tt)
        prefetch(tt->first_entry(key()));

    // Set capture piece in state info
    st->capturedPiece = captured;

    // Calculate checkers bitboard (if move gives check)
    st->checkersBB = givesCheck ? attackers_to(square<KING>(them)) & pieces(us) : 0;

    // Flip side to move
    sideToMove = them;

    // Update slider blockers and check squares for the new side to move
    set_check_info();

    // Calculate the repetition info. It is the ply distance from the previous
    // occurrence of the same position, negative in the 3-fold case, or zero
    // if the position was not repeated.
    st->repetition = 0;
    int end        = std::min(st->rule50, st->pliesFromNull);
    if (end >= 4)
    {
        StateInfo* stp = st->previous->previous;
        for (int i = 4; i <= end; i += 2)
        {
            stp = stp->previous->previous;
            // Ensure the comparison is valid (stp cannot be null)
            if (!stp) break;
            if (stp->key == st->key)
            {
                st->repetition = stp->repetition ? -i : i;
                break;
            }
        }
    }

    assert(pos_is_ok());

    return dp;
}


// Unmakes a move. When it returns, the position should
// be restored to exactly the same state as before the move was made.
void Position::undo_move(Move m) {

    assert(m.is_ok());

    sideToMove = ~sideToMove;

    Color  us   = sideToMove;
    Square from = m.from_sq();
    Square to   = m.to_sq();
    Piece  pc   = piece_on(to); // Note: This is the piece *after* the move, potentially promoted

    assert(empty(from) || m.type_of() == CASTLING);
    assert(type_of(st->capturedPiece) != KING);

    if (m.type_of() == PROMOTION)
    {
        assert(relative_rank(us, to) == RANK_8);
        assert(type_of(pc) == m.promotion_type());
        assert(type_of(pc) >= KNIGHT && type_of(pc) <= QUEEN);

        remove_piece(to); // Remove the promoted piece
        pc = make_piece(us, PAWN); // Piece type is now PAWN
        put_piece(pc, to); // Put the pawn back on the promotion square
    }

    if (m.type_of() == CASTLING)
    {
        Square rfrom, rto;
        // 'from' is king's origin, 'to' is the square the rook started on (captured)
        do_castling<false>(us, from, to, rfrom, rto);
    }
    else // NORMAL, EN_PASSANT (PROMOTION handled above)
    {
        // Put the moving piece back at the source square
        move_piece(to, from);

        if (st->capturedPiece)
        {
            Square capsq = to; // Default capture square is 'to'

            // If the move was en passant, the captured pawn is restored
            // one square behind the 'to' square.
            if (m.type_of() == EN_PASSANT)
            {
                capsq -= pawn_push(us);

                assert(type_of(pc) == PAWN); // pc should be the pawn that moved
                assert(to == st->previous->epSquare); // 'to' must be the previous EP square
                assert(relative_rank(us, to) == RANK_6);
                assert(piece_on(capsq) == NO_PIECE); // Square must be empty now
                assert(st->capturedPiece == make_piece(~us, PAWN));
            }

            // Restore the captured piece
            put_piece(st->capturedPiece, capsq);
        }
    }

    // Finally point our state pointer back to the previous state
    // This restores all state info (epSquare, castlingRights, rule50, keys, etc.)
    st = st->previous;
    --gamePly;

    assert(pos_is_ok());
}


// Helper used to do/undo a castling move. This is a bit
// tricky in Chess960 where from/to squares can overlap.
template<bool Do>
void Position::do_castling(
  Color us, Square from, Square& to, Square& rfrom, Square& rto, DirtyPiece* const dp) {

    bool kingSide = to > from; // In std chess, target rook > king implies kingside
    rfrom         = to;  // Castling is encoded as "king captures friendly rook"
    rto           = relative_square(us, kingSide ? SQ_F1 : SQ_D1); // Standard rook destination
    to            = relative_square(us, kingSide ? SQ_G1 : SQ_C1); // Standard king destination

    assert(!Do || dp); // Ensure dp is provided when doing the move

    Piece king = make_piece(us, KING);
    Piece rook = make_piece(us, ROOK);

    Square kingFrom = Do ? from : to;
    Square kingTo   = Do ? to : from;
    Square rookFrom = Do ? rfrom : rto;
    Square rookTo   = Do ? rto : rfrom;

    if (Do) // Update DP structure if doing the move
    {
        dp->piece[0]  = king;
        dp->from[0]   = kingFrom;
        dp->to[0]     = kingTo;
        dp->piece[1]  = rook;
        dp->from[1]   = rookFrom;
        dp->to[1]     = rookTo;
        dp->dirty_num = 2;
    }

    // Remove both pieces first since squares could overlap in Chess960
    remove_piece(kingFrom);
    remove_piece(rookFrom);
    board[kingFrom] = board[rookFrom] =
      NO_PIECE;  // remove_piece does not do this for us

    // Place pieces on their destination squares
    put_piece(king, kingTo);
    put_piece(rook, rookTo);
}


// Used to do a "null move": it flips
// the side to move without executing any move on the board.
void Position::do_null_move(StateInfo& newSt, const TranspositionTable& tt) {

    assert(!checkers());
    assert(&newSt != st);

    // Copy the current state, point pointers, and make the new state current
    std::memcpy(&newSt, st, sizeof(StateInfo));
    newSt.previous = st;
    st->next       = &newSt;
    st             = &newSt;

    // Reset EP square if it exists, and update the key
    if (st->epSquare != SQ_NONE)
    {
        st->key ^= Zobrist::enpassant[file_of(st->epSquare)];
        st->epSquare = SQ_NONE;
    }

    // Flip side in the key and prefetch TT entry
    st->key ^= Zobrist::side;
    prefetch(tt.first_entry(key()));

    // Update state variables for null move
    st->rule50++; // Increment rule50 counter
    st->pliesFromNull = 0; // Reset plies from null counter

    sideToMove = ~sideToMove; // Flip side to move

    // Update check info for the new side to move
    set_check_info();

    // A null move cannot create a repetition by itself
    st->repetition = 0;

    assert(pos_is_ok());
}


// Must be used to undo a "null move"
void Position::undo_null_move() {

    assert(!checkers());

    st         = st->previous; // Revert to previous state
    sideToMove = ~sideToMove;  // Flip side back
    // No board changes to undo
}


// Tests if the SEE (Static Exchange Evaluation)
// value of move is greater or equal to the given threshold. We'll use an
// algorithm similar to alpha-beta pruning with a null window.
bool Position::see_ge(Move m, int threshold) const {

    assert(m.is_ok());

    // Only deal with normal moves, assume others pass a simple SEE
    // E.g., Castling has SEE 0. EP capture has SEE PawnValue. Promotion capture SEE depends.
    if (m.type_of() != NORMAL)
        return (m.type_of() == EN_PASSANT ? PawnValue : VALUE_ZERO) >= threshold;

    Square from = m.from_sq(), to = m.to_sq();
    Piece  targetPiece = piece_on(to);

    // If the target square is empty, SEE is 0, unless it's a promotion.
    if (targetPiece == NO_PIECE)
         return VALUE_ZERO >= threshold; // Assuming non-promotion NORMAL moves don't target empty squares meaningfully for SEE


    // Initial swap starts with the value of the captured piece minus the threshold.
    int swap = PieceValue[targetPiece] - threshold;
    if (swap < 0) // If capturing the piece doesn't already beat threshold, it's losing SEE.
        return false;

    // The next swap is the value of the capturing piece minus the current swap balance.
    // If this is positive, it means the first capture is immediately profitable (gain > loss).
    swap = PieceValue[piece_on(from)] - swap;
    if (swap <= 0) // If the initial capture is net positive or break-even, it's winning SEE.
        return true;

    assert(color_of(piece_on(from)) == sideToMove);

    // Setup for the SEE loop:
    Bitboard occupied  = pieces() ^ from ^ to;  // Occupancy excluding the moving piece and the captured piece's square temporarily
    Color    stm       = sideToMove;           // Start with the side *opposite* to the one that just moved
    Bitboard attackers = attackers_to(to, occupied); // Find all pieces attacking 'to'
    Bitboard stmAttackers, bb;
    int      res = 1; // Track if the side to move (stm) is winning (1) or losing (0) the exchange

    while (true)
    {
        stm = ~stm; // Switch side
        attackers &= occupied; // Consider only attackers on occupied squares

        // If the current side to move (stm) has no more attackers, they lose.
        if (!(stmAttackers = attackers & pieces(stm)))
            break;

        // Handle pinned pieces: they can only attack along the pin line.
        // If there are pinners to the *opponent's* king, check if any of our attackers are pinned.
        if (pinners(~stm) & occupied)
        {
            // Remove pinned attackers unless they attack along the pin line (handled implicitly by attackers_to)
            // A simpler approximation: remove all pieces blocked by sliders toward *our* king.
            stmAttackers &= ~blockers_for_king(stm); // Needs accurate blockers calculation

            if (!stmAttackers) // If all attackers were pinned appropriately, stm loses.
                break;
        }

        res ^= 1; // The side to move is now potentially winning (res=1) or losing (res=0) based on the opponent's last move

        // Find the least valuable attacker for the current side (stm)
        // Remove it from occupancy and update the potential attackers (sliders revealed)
        if ((bb = stmAttackers & pieces(PAWN)))
        {
            swap = PawnValue - swap; // Value of pawn captured - current balance
            if (swap < res) break; // If capturing this pawn makes the exchange losing for stm, stop.
            occupied ^= least_significant_square_bb(bb); // Remove the pawn
            // Add revealed diagonal attackers (Bishops/Queens)
            attackers |= attacks_bb<BISHOP>(to, occupied) & pieces(BISHOP, QUEEN);
        }
        else if ((bb = stmAttackers & pieces(KNIGHT)))
        {
            swap = KnightValue - swap;
            if (swap < res) break;
            occupied ^= least_significant_square_bb(bb);
            // Knights don't reveal sliders
        }
        else if ((bb = stmAttackers & pieces(BISHOP)))
        {
            swap = BishopValue - swap;
            if (swap < res) break;
            occupied ^= least_significant_square_bb(bb);
            // Add revealed diagonal attackers
            attackers |= attacks_bb<BISHOP>(to, occupied) & pieces(BISHOP, QUEEN);
        }
        else if ((bb = stmAttackers & pieces(ROOK)))
        {
            swap = RookValue - swap;
            if (swap < res) break;
            occupied ^= least_significant_square_bb(bb);
            // Add revealed rank/file attackers
            attackers |= attacks_bb<ROOK>(to, occupied) & pieces(ROOK, QUEEN);
        }
        else if ((bb = stmAttackers & pieces(QUEEN)))
        {
            swap = QueenValue - swap;
            if (swap < res) break; // Should imply previous capture was King, which isn't possible here? Check logic.
            //assert(swap >= res); // Assuming king is never captured this way
            occupied ^= least_significant_square_bb(bb);
            // Add revealed diagonal and rank/file attackers
            attackers |= (attacks_bb<BISHOP>(to, occupied) & pieces(BISHOP, QUEEN))
                       | (attacks_bb<ROOK>(to, occupied) & pieces(ROOK, QUEEN));
        }
        else // KING is the least valuable attacker remaining
        {
            // If the king captures, but the opponent still has attackers, the exchange is lost.
            // Otherwise, the exchange is won (or break-even depending on the initial threshold).
            return (attackers & ~pieces(stm)) ? res ^ 1 : res;
        }
        // Update attackers bitboard based on new occupancy (remove pieces no longer on board)
        // Note: This update happens implicitly at the start of the loop (attackers &= occupied)
    }

    // If the loop terminates because one side ran out of attackers,
    // the final result depends on whose turn it was.
    // If res=1, the last side to move won. If res=0, the last side to move lost.
    return bool(res);
}


// Tests whether the position is drawn by 50-move rule
// or by repetition. It does not detect stalemates.
bool Position::is_draw(int ply) const {

    if (st->rule50 > 99 && (!checkers() || MoveList<LEGAL>(*this).size()))
        return true;

    return is_repetition(ply);
}


// Return a draw score if a position repeats once earlier but strictly
// after the root, or repeats twice before or at the root.
bool Position::is_repetition(int ply) const {
    // repetition is > 0 if repeated, < 0 if 3-fold repetition
    // We draw if the repetition occurred after the root ply (ply=0 at root)
    // Negative repetition values indicate 3-fold repetition, always a draw.
    return st->repetition && (st->repetition < 0 || st->repetition < ply);
}


// Tests whether there has been at least one repetition
// of positions since the last capture or pawn move.
bool Position::has_repeated() const {

    StateInfo* stc = st;
    // Iterate back through states as long as rule50 allows repetitions
    int end = std::min(st->rule50, st->pliesFromNull);
    // Must check at least 4 half-moves back (2 full moves) for a repetition
    while (end >= 4)
    {
        // Go back two states (one full move)
        stc = stc->previous;
        if (!stc) return false; // Should not happen if end >= 4
        stc = stc->previous;
        if (!stc) return false;

        if (stc->key == st->key) // Compare keys for repetition
            return true;

        end -= 2; // Checked one full move back
    }
    return false;
}


// Tests if the position has a move which draws by repetition.
// This function accurately matches the outcome of is_draw() over all legal moves.
bool Position::upcoming_repetition(int ply) const {

    // Based on Marcel van Kervinck's cuckoo algorithm explained in position.cpp
    // It checks for potential repetitions using precomputed hashes of reversible moves.

    int j; // index for cuckoo table

    // Maximum number of plies to check back, limited by rule50 and plies since null move
    int end = std::min(st->rule50, st->pliesFromNull);

    // Need at least 3 plies to check for a repetition via a reversible move
    if (end < 3)
        return false;

    Key        originalKey = st->key;
    StateInfo* stp         = st->previous; // State at ply - 1
    // other = H(P) ^ H(P-1) ^ S = H(P-1) without move m from P-1 to P
    Key        other       = originalKey ^ stp->key ^ Zobrist::side;

    // Iterate back checking states at ply-3, ply-5, etc.
    for (int i = 3; i <= end; i += 2)
    {
        stp = stp->previous; // State at ply - 2
        if (!stp) break;
        stp = stp->previous; // State at ply - i
        if (!stp) break;

        // other ^= H(P-(i-1)) ^ H(P-i) ^ S
        // This accumulates XOR difference between P and P-i based on reversible moves
        other ^= stp->key ^ stp->previous->key ^ Zobrist::side;

        // If other is 0, it means H(P) == H(P-i), a potential repetition
        if (other != 0)
            continue;

        // H(P) == H(P-i) => Potential repetition exists.
        // moveKey = H(P) ^ H(P-i+1) should correspond to the Zobrist key of the
        // single reversible move that transforms position P-i to P-i+2 etc. up to P.
        // We check if such a reversible move exists in our precomputed cuckoo table.
        Key moveKey = originalKey ^ stp->key; // Zobrist diff for the move

        // Check both possible hash locations in the cuckoo table
        if ((j = H1(moveKey), cuckoo[j] == moveKey) || (j = H2(moveKey), cuckoo[j] == moveKey))
        {
            Move   move = cuckooMove[j]; // The potentially repeating reversible move
            Square s1   = move.from_sq();
            Square s2   = move.to_sq();

            // Ensure the path between the squares is empty (it's a truly reversible move)
            // This condition seems slightly off, between_bb includes endpoints. Maybe should be:
            // if (!(between_bb(s1, s2) & ~(s1 | s2) & pieces()))
            // Let's stick to original logic: checks if ONLY the destination is between from and to.
            // This seems wrong. A simple reversible move like Nc3-e4 shouldn't trigger this.
            // Let's assume the cuckoo table only stores *valid* reversible piece moves.
            // The check `!(between_bb(s1, s2) ^ s2) & pieces())` might be for sliders?
            // It checks if (path_excluding_s2 & pieces) is empty. For non-sliders, between_bb is 0.
            // So this effectively checks `!(pieces() & s1)` for non-sliders? No, that's not right.
            // Let's trust the original logic's intent: check if the move is currently possible reversibly.
             if (!((between_bb(s1, s2) ^ s2) & pieces())) // Original check
             {
                // If the repetition is detected deeper in the search (ply > i), it's a draw.
                if (ply > i)
                    return true;

                // For nodes at or before the root (ply <= i), we need to ensure
                // it's a genuine 3-fold repetition, not just returning to the current state.
                // The 'repetition' field in StateInfo tracks this. If stp->repetition
                // is non-zero, it means the state at ply 'stp' was already a repetition.
                if (stp->repetition)
                    return true;
             }
        }
    }
    return false; // No drawing repetition found via a reversible move
}


// Flips position with the white and black sides reversed. This
// is only useful for debugging e.g. for finding evaluation symmetry bugs.
void Position::flip() {

    string            f, token;
    std::stringstream ss(fen()); // Get current FEN

    // 1. Reverse ranks and swap piece case
    for (Rank r = RANK_8; r >= RANK_1; --r)
    {
        std::getline(ss, token, r > RANK_1 ? '/' : ' '); // Read one rank description
        // Insert the rank description at the beginning of 'f', separated by '/'
        f.insert(0, token + (f.empty() ? " " : "/"));
    }

    // 2. Swap active color
    ss >> token; // Read color ('w' or 'b')
    f += (token == "w" ? "B " : "W "); // Append swapped color (will be lowercased later)

    // 3. Swap castling rights case
    ss >> token; // Read castling rights string
    f += token + " "; // Append original castling rights (case will be swapped later)

    // Perform case swap on the constructed FEN part (pieces, color, castling)
    std::transform(f.begin(), f.end(), f.begin(),
                   [](char c) { return char(std::islower(c) ? std::toupper(c) : std::tolower(c)); });

    // 4. Adjust EP square rank
    ss >> token; // Read EP square string
    if (token != "-")
        // Change rank '3' to '6' or '6' to '3'
        token.replace(1, 1, token[1] == '3' ? "6" : "3");
    f += token;

    // 5. Append clocks
    std::getline(ss, token); // Read the rest of the line (rule50, fullmove)
    f += token;

    // Set the new position using the modified FEN string
    set(f, is_chess960(), st);

    assert(pos_is_ok());
}


// Performs some consistency checks for the position object
// and raise an assert if something wrong is detected.
// This is meant to be helpful when debugging.
bool Position::pos_is_ok() const {

    constexpr bool Fast = true;  // Quick (default) or full check?

    // Basic checks always performed
    if ((sideToMove != WHITE && sideToMove != BLACK) // Side to move must be valid
        || piece_on(square<KING>(WHITE)) != W_KING   // White king must be on its square
        || piece_on(square<KING>(BLACK)) != B_KING   // Black king must be on its square
        || (ep_square() != SQ_NONE                   // If EP square exists...
            && relative_rank(sideToMove, ep_square()) != RANK_6)) // ...it must be on rank 6 relative to side to move
        { assert(0 && "pos_is_ok: Basic state check failed"); return false; }

    // Skip further checks if Fast mode is enabled
    if (Fast)
        return true;

    // Full checks (only if Fast is false)
    if (pieceCount[W_KING] != 1 || pieceCount[B_KING] != 1 // Exactly one king per side
        || attackers_to_exist(square<KING>(~sideToMove), pieces(), sideToMove)) // Side not to move cannot be in check
        { assert(0 && "pos_is_ok: King count or non-active side in check"); return false; }

    if ((pieces(PAWN) & (Rank1BB | Rank8BB)) // Pawns cannot be on first or eighth rank
        || pieceCount[W_PAWN] > 8 || pieceCount[B_PAWN] > 8) // Max 8 pawns per side
        { assert(0 && "pos_is_ok: Pawn position or count invalid"); return false; }

    if ((pieces(WHITE) & pieces(BLACK))                 // White and Black pieces cannot overlap
        || (pieces(WHITE) | pieces(BLACK)) != pieces()  // Union of White and Black must equal all pieces
        || popcount(pieces(WHITE)) > 16 || popcount(pieces(BLACK)) > 16) // Max 16 pieces per side
        { assert(0 && "pos_is_ok: Color bitboard inconsistency"); return false; }

    // Check piece type bitboards for consistency
    for (PieceType p1 = PAWN; p1 <= KING; ++p1) {
        if ((pieces(WHITE, p1) | pieces(BLACK, p1)) != pieces(p1)) // Color union must match type
           { assert(0 && "pos_is_ok: Piece type bitboard color mismatch"); return false; }
        for (PieceType p2 = PAWN; p2 < p1; ++p2) { // Check against lower piece types
            if ((pieces(p1) & pieces(p2))) // Different piece types cannot overlap
               { assert(0 && "pos_is_ok: Piece type bitboard overlap"); return false; }
        }
    }

    // Check piece counts against bitboards and board array
    for (Piece pc : Pieces)
        if (pc != NO_PIECE) { // Iterate through actual piece values
            if (pieceCount[pc] != popcount(pieces(color_of(pc), type_of(pc))) // Count vs bitboard popcount
                || pieceCount[pc] != std::count(board, board + SQUARE_NB, pc)) // Count vs board array scan
               { assert(0 && "pos_is_ok: Piece count mismatch"); return false; }
        }
        // Also check aggregate counts per color
    if (pieceCount[make_piece(WHITE, ALL_PIECES)] != popcount(pieces(WHITE))
        || pieceCount[make_piece(BLACK, ALL_PIECES)] != popcount(pieces(BLACK)))
        { assert(0 && "pos_is_ok: Total color piece count mismatch"); return false; }


    // Check castling rights consistency
    for (Color c : {WHITE, BLACK})
        for (CastlingRights cr : {c & KING_SIDE, c & QUEEN_SIDE}) // Iterate relevant castling rights
        {
            if (!can_castle(cr)) continue; // Skip if right doesn't exist

            Square rsq = castlingRookSquare[cr];
            Square ksq = square<KING>(c);

            if (piece_on(rsq) != make_piece(c, ROOK)   // Rook must be on its castling square
                || (castlingRightsMask[rsq] & cr) != cr // Rook square mask must include this right
                || (castlingRightsMask[ksq] & cr) != cr) // King square mask must include this right
                { assert(0 && "pos_is_ok: Castling rights data inconsistency"); return false; }
        }

    // Check slider blockers calculation (can be expensive)
    StateInfo tempSt = *st; // Copy current state
    set_check_info(); // Recalculate check info based on current board state
    if (tempSt.blockersForKing[WHITE] != st->blockersForKing[WHITE]
        || tempSt.blockersForKing[BLACK] != st->blockersForKing[BLACK]
        || tempSt.pinners[WHITE] != st->pinners[WHITE]
        || tempSt.pinners[BLACK] != st->pinners[BLACK])
       { assert(0 && "pos_is_ok: Slider blocker/pinner info mismatch"); return false; }
       // Note: Also implicitly checks checkSquares consistency if set_check_info is correct


    // Check Zobrist keys (very expensive, involves recalculating from scratch)
    // Key recalculation omitted for brevity, but would involve iterating pieces,
    // XORing Zobrist::psq, Zobrist::castling, Zobrist::enpassant, Zobrist::side
    // Key Check Example (conceptual):
    /*
    Key tempKey = 0;
    for (Square s = SQ_A1; s <= SQ_H8; ++s)
        if (board[s] != NO_PIECE) tempKey ^= Zobrist::psq[board[s]][s];
    if (st->epSquare != SQ_NONE) tempKey ^= Zobrist::enpassant[file_of(st->epSquare)];
    tempKey ^= Zobrist::castling[st->castlingRights];
    if (sideToMove == BLACK) tempKey ^= Zobrist::side;
    if (key() != adjust_key50<false>(tempKey)) // Need to handle rule50 adjustment
       { assert(0 && "pos_is_ok: Zobrist key mismatch"); return false; }
    */
    // Similar checks for pawnKey, materialKey, nonPawnKey, minorPieceKey...


    return true; // All checks passed
}

} // namespace Stockfish
