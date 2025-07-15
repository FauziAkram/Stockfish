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

inline int H1(Key h) { return h & 0x1fff; }
inline int H2(Key h) { return (h >> 16) & 0x1fff; }
std::array<Key, 8192>  cuckoo;
std::array<Move, 8192> cuckooMove;

void Position::init() {
    PRNG rng(1070372);
    for (Piece pc : Pieces)
        for (Square s = SQ_A1; s <= SQ_H8; ++s)
            Zobrist::psq[pc][s] = rng.rand<Key>();
    for (File f = FILE_A; f <= FILE_H; ++f)
        Zobrist::enpassant[f] = rng.rand<Key>();
    for (int cr = NO_CASTLING; cr <= ANY_CASTLING; ++cr)
        Zobrist::castling[cr] = rng.rand<Key>();
    Zobrist::side    = rng.rand<Key>();
    Zobrist::noPawns = rng.rand<Key>();
    cuckoo.fill(0);
    cuckooMove.fill(Move::none());
    [[maybe_unused]] int count = 0;
    for (Piece pc : Pieces)
        for (Square s1 = SQ_A1; s1 <= SQ_H8; ++s1)
            for (Square s2 = Square(s1 + 1); s2 <= SQ_H8; ++s2)
                if ((type_of(pc) != PAWN) && (attacks_bb<KING>(s1) & s2)) // Simplified check
                {
                    Move move = Move(s1, s2);
                    Key  key  = Zobrist::psq[pc][s1] ^ Zobrist::psq[pc][s2] ^ Zobrist::side;
                    int  i    = H1(key);
                    while (true)
                    {
                        std::swap(cuckoo[i], key);
                        std::swap(cuckooMove[i], move);
                        if (move == Move::none())
                            break;
                        i = (i == H1(key)) ? H2(key) : H1(key);
                    }
                    count++;
                }
}

Position& Position::set(const string& fenStr, bool isChess960, StateInfo* si) {
    unsigned char      col, row, token;
    size_t             idx;
    Square             sq = SQ_A8;
    std::istringstream ss(fenStr);
    std::memset(this, 0, sizeof(Position));
    std::memset(si, 0, sizeof(StateInfo));
    st = si;
    ss >> std::noskipws;
    while ((ss >> token) && !isspace(token))
    {
        if (isdigit(token))
            sq += (token - '0') * EAST;
        else if (token == '/')
            sq += 2 * SOUTH;
        else if ((idx = PieceToChar.find(token)) != string::npos)
        {
            put_piece(Piece(idx), sq);
            ++sq;
        }
    }
    ss >> token;
    sideToMove = (token == 'w' ? WHITE : BLACK);
    ss >> token;
    while ((ss >> token) && !isspace(token))
    {
        Square rsq;
        Color  c    = islower(token) ? BLACK : WHITE;
        Piece  rook = make_piece(c, ROOK);
        token = char(toupper(token));
        if (token == 'K')
            for (rsq = relative_square(c, SQ_H1); piece_on(rsq) != rook; --rsq) {}
        else if (token == 'Q')
            for (rsq = relative_square(c, SQ_A1); piece_on(rsq) != rook; ++rsq) {}
        else if (token >= 'A' && token <= 'H')
            rsq = make_square(File(token - 'A'), relative_rank(c, RANK_1));
        else
            continue;
        set_castling_right(c, rsq);
    }
    bool enpassant = false;
    if (((ss >> col) && (col >= 'a' && col <= 'h'))
        && ((ss >> row) && (row == (sideToMove == WHITE ? '6' : '3'))))
    {
        st->epSquare = make_square(File(col - 'a'), Rank(row - '1'));
        enpassant = (attacks_bb<PAWN>(st->epSquare, ~sideToMove) & pieces(sideToMove, PAWN))
                 && (pieces(~sideToMove, PAWN) & (st->epSquare + pawn_push(~sideToMove)))
                 && !(pieces() & (st->epSquare | (st->epSquare + pawn_push(sideToMove))));
    }
    if (!enpassant)
        st->epSquare = SQ_NONE;
    ss >> std::skipws >> st->rule50 >> gamePly;
    gamePly = std::max(2 * (gamePly - 1), 0) + (sideToMove == BLACK);
    chess960 = isChess960;
    set_state();
    assert(pos_is_ok());
    return *this;
}

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

void Position::set_state() const {
    st->key = st->materialKey = 0;
    st->pawnKey = Zobrist::noPawns;
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
        else if (type_of(pc) != KING)
            st->nonPawnMaterial[color_of(pc)] += PieceValue[pc];
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

Position& Position::set(const string& code, Color c, StateInfo* si) {
    assert(code[0] == 'K');
    string sides[] = {code.substr(code.find('K', 1)),
                      code.substr(0, std::min(code.find('v'), code.find('K', 1)))};
    assert(sides[0].length() > 0 && sides[0].length() < 8);
    assert(sides[1].length() > 0 && sides[1].length() < 8);
    std::transform(sides[c].begin(), sides[c].end(), sides[c].begin(), ::tolower);
    string fenStr = "8/" + sides[0] + char(8 - sides[0].length() + '0') + "/8/8/8/8/" + sides[1]
                  + char(8 - sides[1].length() + '0') + "/8 w - - 0 10";
    return set(fenStr, false, si);
}

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

void Position::update_slider_blockers(Color c) const {
    Square ksq = square<KING>(c);
    st->blockersForKing[c] = 0;
    st->pinners[~c]        = 0;
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

bool Position::legal(Move m) const {
    assert(m.is_ok());
    Color  us   = sideToMove;
    Square from = m.from_sq();
    Square to   = m.to_sq();
    assert(color_of(moved_piece(m)) == us);
    assert(piece_on(square<KING>(us)) == make_piece(us, KING));
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
    if (m.type_of() == CASTLING)
    {
        to             = relative_square(us, to > from ? SQ_G1 : SQ_C1);
        Direction step = to > from ? WEST : EAST;
        for (Square s = to; s != from; s += step)
            if (attackers_to_exist(s, pieces(), ~us))
                return false;
        return !chess960 || !(blockers_for_king(us) & m.to_sq());
    }
    if (type_of(piece_on(from)) == KING)
        return !(attackers_to_exist(to, pieces() ^ from, ~us));
    return !(blockers_for_king(us) & from) || line_bb(from, to) & pieces(us, KING);
}

bool Position::pseudo_legal(const Move m) const {
    Color  us   = sideToMove;
    Square from = m.from_sq();
    Square to   = m.to_sq();
    Piece  pc   = moved_piece(m);
    if (m.type_of() != NORMAL)
        return checkers() ? MoveList<EVASIONS>(*this).contains(m)
                          : MoveList<NON_EVASIONS>(*this).contains(m);
    assert(m.promotion_type() - KNIGHT == NO_PIECE_TYPE);
    if (pc == NO_PIECE || color_of(pc) != us)
        return false;
    if (pieces(us) & to)
        return false;
    if (type_of(pc) == PAWN)
    {
        if ((Rank8BB | Rank1BB) & to)
            return false;
        const bool isCapture    = bool(attacks_bb<PAWN>(from, us) & pieces(~us) & to);
        const bool isSinglePush = (from + pawn_push(us) == to) && empty(to);
        const bool isDoublePush = (from + 2 * pawn_push(us) == to)
                               && (relative_rank(us, from) == RANK_2) && empty(to)
                               && empty(to - pawn_push(us));
        if (!(isCapture || isSinglePush || isDoublePush))
            return false;
    }
    else if (!(attacks_bb(type_of(pc), from, pieces()) & to))
        return false;
    if (checkers())
    {
        if (type_of(pc) != KING)
        {
            if (more_than_one(checkers()))
                return false;
            if (!(between_bb(square<KING>(us), lsb(checkers())) & to))
                return false;
        }
        else if (attackers_to_exist(to, pieces() ^ from, ~us))
            return false;
    }
    return true;
}

bool Position::gives_check(Move m) const {
    assert(m.is_ok());
    assert(color_of(moved_piece(m)) == sideToMove);
    Square from = m.from_sq();
    Square to   = m.to_sq();
    if (check_squares(type_of(piece_on(from))) & to)
        return true;
    if (blockers_for_king(~sideToMove) & from)
        return !(line_bb(from, to) & pieces(~sideToMove, KING)) || m.type_of() == CASTLING;
    switch (m.type_of())
    {
    case NORMAL: return false;
    case PROMOTION:
        return attacks_bb(m.promotion_type(), to, pieces() ^ from) & pieces(~sideToMove, KING);
    case EN_PASSANT: {
        Square   capsq = make_square(file_of(to), rank_of(from));
        Bitboard b     = (pieces() ^ from ^ capsq) | to;
        return (attacks_bb<ROOK>(square<KING>(~sideToMove), b) & pieces(sideToMove, QUEEN, ROOK))
             | (attacks_bb<BISHOP>(square<KING>(~sideToMove), b) & pieces(sideToMove, QUEEN, BISHOP));
    }
    default: { //CASTLING
        Square rto = relative_square(sideToMove, to > from ? SQ_F1 : SQ_D1);
        return check_squares(ROOK) & rto;
    }}
}

void Position::do_move(Move m, StateInfo& newSt, const TranspositionTable* tt) {
    assert(m.is_ok());
    assert(&newSt != st);
    Key k = st->key ^ Zobrist::side;
    std::memcpy(&newSt, st, offsetof(StateInfo, key));
    newSt.previous = st;
    st             = &newSt;
    ++gamePly;
    ++st->rule50;
    ++st->pliesFromNull;
    Color  us       = sideToMove;
    Color  them     = ~us;
    Square from     = m.from_sq();
    Square to       = m.to_sq();
    Piece  pc       = piece_on(from);
    Piece  captured = m.type_of() == EN_PASSANT ? make_piece(them, PAWN) : piece_on(to);
    assert(color_of(pc) == us);
    assert(captured == NO_PIECE || color_of(captured) == (m.type_of() != CASTLING ? them : us));
    assert(type_of(captured) != KING);
    if (m.type_of() == CASTLING)
    {
        assert(pc == make_piece(us, KING));
        assert(captured == make_piece(us, ROOK));
        Square rfrom, rto;
        do_castling<true>(us, from, to, rfrom, rto);
        k ^= Zobrist::psq[captured][rfrom] ^ Zobrist::psq[captured][rto];
        captured = NO_PIECE;
    }
    if (captured)
    {
        Square capsq = to;
        if (type_of(captured) == PAWN)
        {
            if (m.type_of() == EN_PASSANT)
            {
                capsq -= pawn_push(us);
                assert(pc == make_piece(us, PAWN));
                assert(to == st->epSquare);
                assert(relative_rank(us, to) == RANK_6);
                assert(piece_on(to) == NO_PIECE);
                assert(piece_on(capsq) == make_piece(them, PAWN));
            }
            st->pawnKey ^= Zobrist::psq[captured][capsq];
        }
        else
            st->nonPawnMaterial[them] -= PieceValue[captured];
        remove_piece(capsq);
        k ^= Zobrist::psq[captured][capsq];
        st->materialKey ^= Zobrist::psq[captured][8 + pieceCount[captured]];
        if (tt)
            prefetch(tt->first_entry(st->materialKey));
        st->rule50 = 0;
    }
    k ^= Zobrist::psq[pc][from] ^ Zobrist::psq[pc][to];
    if (st->epSquare != SQ_NONE)
    {
        k ^= Zobrist::enpassant[file_of(st->epSquare)];
        st->epSquare = SQ_NONE;
    }
    if (st->castlingRights && (castlingRightsMask[from] | castlingRightsMask[to]))
    {
        k ^= Zobrist::castling[st->castlingRights];
        st->castlingRights &= ~(castlingRightsMask[from] | castlingRightsMask[to]);
        k ^= Zobrist::castling[st->castlingRights];
    }
    if (m.type_of() != CASTLING)
        move_piece(from, to);
    if (type_of(pc) == PAWN)
    {
        if ((int(to) ^ int(from)) == 16
            && (attacks_bb<PAWN>(to - pawn_push(us), us) & pieces(them, PAWN)))
        {
            st->epSquare = to - pawn_push(us);
            k ^= Zobrist::enpassant[file_of(st->epSquare)];
        }
        else if (m.type_of() == PROMOTION)
        {
            Piece promotion = make_piece(us, m.promotion_type());
            assert(relative_rank(us, to) == RANK_8);
            assert(type_of(promotion) >= KNIGHT && type_of(promotion) <= QUEEN);
            remove_piece(to);
            put_piece(promotion, to);
            k ^= Zobrist::psq[pc][to] ^ Zobrist::psq[promotion][to];
            st->pawnKey ^= Zobrist::psq[pc][to];
            st->materialKey ^= Zobrist::psq[promotion][8 + pieceCount[promotion] - 1]
                             ^ Zobrist::psq[pc][8 + pieceCount[pc]];
            st->nonPawnMaterial[us] += PieceValue[promotion];
        }
        st->pawnKey ^= Zobrist::psq[pc][from] ^ Zobrist::psq[pc][to];
        st->rule50 = 0;
    }
    st->capturedPiece = captured;
    st->key = k;
    st->checkersBB = gives_check(m) ? attackers_to(square<KING>(them)) & pieces(us) : 0;
    sideToMove = ~sideToMove;
    set_check_info();
    st->repetition = 0;
    int end = std::min(st->rule50, st->pliesFromNull);
    if (end >= 4)
    {
        StateInfo* stp = st->previous->previous;
        for (int i = 4; i <= end; i += 2)
        {
            stp = stp->previous->previous;
            if (stp->key == st->key)
            {
                st->repetition = stp->repetition ? -i : i;
                break;
            }
        }
    }
    assert(pos_is_ok());
}

void Position::undo_move(Move m) {
    assert(m.is_ok());
    sideToMove = ~sideToMove;
    Color  us   = sideToMove;
    Square from = m.from_sq();
    Square to   = m.to_sq();
    Piece  pc   = piece_on(to);
    assert(empty(from) || m.type_of() == CASTLING);
    assert(type_of(st->capturedPiece) != KING);
    if (m.type_of() == PROMOTION)
    {
        assert(relative_rank(us, to) == RANK_8);
        assert(type_of(pc) == m.promotion_type());
        assert(type_of(pc) >= KNIGHT && type_of(pc) <= QUEEN);
        remove_piece(to);
        pc = make_piece(us, PAWN);
        put_piece(pc, to);
    }
    if (m.type_of() == CASTLING)
    {
        Square rfrom, rto;
        do_castling<false>(us, from, to, rfrom, rto);
    }
    else
    {
        move_piece(to, from);
        if (st->capturedPiece)
        {
            Square capsq = to;
            if (m.type_of() == EN_PASSANT)
            {
                capsq -= pawn_push(us);
                assert(type_of(pc) == PAWN);
                assert(to == st->previous->epSquare);
                assert(relative_rank(us, to) == RANK_6);
                assert(piece_on(capsq) == NO_PIECE);
                assert(st->capturedPiece == make_piece(~us, PAWN));
            }
            put_piece(st->capturedPiece, capsq);
        }
    }
    st = st->previous;
    --gamePly;
    assert(pos_is_ok());
}

template<bool Do>
void Position::do_castling(Color us, Square from, Square& to, Square& rfrom, Square& rto) {
    bool kingSide = to > from;
    rfrom         = to;
    rto           = relative_square(us, kingSide ? SQ_F1 : SQ_D1);
    to            = relative_square(us, kingSide ? SQ_G1 : SQ_C1);
    remove_piece(Do ? from : to);
    remove_piece(Do ? rfrom : rto);
    board[Do ? from : to] = board[Do ? rfrom : rto] = NO_PIECE;
    put_piece(make_piece(us, KING), Do ? to : from);
    put_piece(make_piece(us, ROOK), Do ? rto : rfrom);
}

void Position::do_null_move(StateInfo& newSt, const TranspositionTable& tt) {
    assert(!checkers());
    assert(&newSt != st);
    std::memcpy(&newSt, st, sizeof(StateInfo));
    newSt.previous = st;
    st             = &newSt;
    if (st->epSquare != SQ_NONE)
    {
        st->key ^= Zobrist::enpassant[file_of(st->epSquare)];
        st->epSquare = SQ_NONE;
    }
    st->key ^= Zobrist::side;
    prefetch(tt.first_entry(key()));
    st->pliesFromNull = 0;
    sideToMove = ~sideToMove;
    set_check_info();
    st->repetition = 0;
    assert(pos_is_ok());
}

void Position::undo_null_move() {
    assert(!checkers());
    st         = st->previous;
    sideToMove = ~sideToMove;
}

bool Position::see_ge(Move m, int threshold) const {
    assert(m.is_ok());
    if (m.type_of() != NORMAL)
        return VALUE_ZERO >= threshold;
    Square from = m.from_sq(), to = m.to_sq();
    int swap = PieceValue[piece_on(to)] - threshold;
    if (swap < 0)
        return false;
    swap = PieceValue[piece_on(from)] - swap;
    if (swap <= 0)
        return true;
    assert(color_of(piece_on(from)) == sideToMove);
    Bitboard occupied  = pieces() ^ from ^ to;
    Color    stm       = sideToMove;
    Bitboard attackers = attackers_to(to, occupied);
    Bitboard stmAttackers, bb;
    int      res = 1;
    while (true)
    {
        stm = ~stm;
        attackers &= occupied;
        if (!(stmAttackers = attackers & pieces(stm)))
            break;
        if (pinners(~stm) & occupied)
        {
            stmAttackers &= ~blockers_for_king(stm);
            if (!stmAttackers)
                break;
        }
        res ^= 1;
        if ((bb = stmAttackers & pieces(PAWN)))
        {
            if ((swap = PawnValue - swap) < res)
                break;
            occupied ^= least_significant_square_bb(bb);
            attackers |= attacks_bb<BISHOP>(to, occupied) & pieces(BISHOP, QUEEN);
        }
        else if ((bb = stmAttackers & pieces(KNIGHT)))
        {
            if ((swap = KnightValue - swap) < res)
                break;
            occupied ^= least_significant_square_bb(bb);
        }
        else if ((bb = stmAttackers & pieces(BISHOP)))
        {
            if ((swap = BishopValue - swap) < res)
                break;
            occupied ^= least_significant_square_bb(bb);
            attackers |= attacks_bb<BISHOP>(to, occupied) & pieces(BISHOP, QUEEN);
        }
        else if ((bb = stmAttackers & pieces(ROOK)))
        {
            if ((swap = RookValue - swap) < res)
                break;
            occupied ^= least_significant_square_bb(bb);
            attackers |= attacks_bb<ROOK>(to, occupied) & pieces(ROOK, QUEEN);
        }
        else if ((bb = stmAttackers & pieces(QUEEN)))
        {
            swap = QueenValue - swap;
            assert(swap >= res);
            occupied ^= least_significant_square_bb(bb);
            attackers |= (attacks_bb<BISHOP>(to, occupied) & pieces(BISHOP, QUEEN))
                       | (attacks_bb<ROOK>(to, occupied) & pieces(ROOK, QUEEN));
        }
        else
            return (attackers & ~pieces(stm)) ? res ^ 1 : res;
    }
    return bool(res);
}

bool Position::is_draw(int ply) const {
    if (st->rule50 > 99 && (!checkers() || MoveList<LEGAL>(*this).size()))
        return true;
    return st->repetition && st->repetition < ply;
}

bool Position::has_repeated() const {
    StateInfo* stc = st;
    int        end = std::min(st->rule50, st->pliesFromNull);
    while (end-- >= 4)
    {
        if (stc->repetition)
            return true;
        stc = stc->previous;
    }
    return false;
}

bool Position::upcoming_repetition(int ply) const {
    int j;
    int end = std::min(st->rule50, st->pliesFromNull);
    if (end < 3)
        return false;
    Key        originalKey = st->key;
    StateInfo* stp         = st->previous;
    Key        other       = originalKey ^ stp->key ^ Zobrist::side;
    for (int i = 3; i <= end; i += 2)
    {
        stp = stp->previous;
        other ^= stp->key ^ stp->previous->key ^ Zobrist::side;
        stp = stp->previous;
        if (other != 0)
            continue;
        Key moveKey = originalKey ^ stp->key;
        if ((j = H1(moveKey), cuckoo[j] == moveKey) || (j = H2(moveKey), cuckoo[j] == moveKey))
        {
            Move   move = cuckooMove[j];
            Square s1   = move.from_sq();
            Square s2   = move.to_sq();
            if (!((between_bb(s1, s2) ^ s2) & pieces()))
            {
                if (ply > i)
                    return true;
                if (stp->repetition)
                    return true;
            }
        }
    }
    return false;
}

void Position::flip() {
    string            f, token;
    std::stringstream ss(fen());
    for (Rank r = RANK_8; r >= RANK_1; --r)
    {
        std::getline(ss, token, r > RANK_1 ? '/' : ' ');
        f.insert(0, token + (f.empty() ? " " : "/"));
    }
    ss >> token;
    f += (token == "w" ? "B " : "W ");
    ss >> token;
    f += token + " ";
    std::transform(f.begin(), f.end(), f.begin(),
                   [](char c) { return char(islower(c) ? toupper(c) : tolower(c)); });
    ss >> token;
    f += (token == "-" ? token : token.replace(1, 1, token[1] == '3' ? "6" : "3"));
    std::getline(ss, token);
    f += token;
    set(f, is_chess960(), st);
    assert(pos_is_ok());
}

bool Position::pos_is_ok() const {
    constexpr bool Fast = true;
    if ((sideToMove != WHITE && sideToMove != BLACK) || piece_on(square<KING>(WHITE)) != W_KING
        || piece_on(square<KING>(BLACK)) != B_KING
        || (ep_square() != SQ_NONE && relative_rank(sideToMove, ep_square()) != RANK_6))
        assert(0 && "pos_is_ok: Default");
    if (Fast)
        return true;
    if (pieceCount[W_KING] != 1 || pieceCount[B_KING] != 1
        || attackers_to_exist(square<KING>(~sideToMove), pieces(), sideToMove))
        assert(0 && "pos_is_ok: Kings");
    if ((pieces(PAWN) & (Rank1BB | Rank8BB)) || pieceCount[W_PAWN] > 8 || pieceCount[B_PAWN] > 8)
        assert(0 && "pos_is_ok: Pawns");
    if ((pieces(WHITE) & pieces(BLACK)) || (pieces(WHITE) | pieces(BLACK)) != pieces()
        || popcount(pieces(WHITE)) > 16 || popcount(pieces(BLACK)) > 16)
        assert(0 && "pos_is_ok: Bitboards");
    for (PieceType p1 = PAWN; p1 <= KING; ++p1)
        for (PieceType p2 = PAWN; p2 <= KING; ++p2)
            if (p1 != p2 && (pieces(p1) & pieces(p2)))
                assert(0 && "pos_is_ok: Bitboards");
    for (Piece pc : Pieces)
        if (pieceCount[pc] != popcount(pieces(color_of(pc), type_of(pc)))
            || pieceCount[pc] != std::count(board, board + SQUARE_NB, pc))
            assert(0 && "pos_is_ok: Pieces");
    for (Color c : {WHITE, BLACK})
        for (CastlingRights cr : {c & KING_SIDE, c & QUEEN_SIDE})
        {
            if (!can_castle(cr))
                continue;
            if (piece_on(castlingRookSquare[cr]) != make_piece(c, ROOK)
                || castlingRightsMask[castlingRookSquare[cr]] != cr
                || (castlingRightsMask[square<KING>(c)] & cr) != cr)
                assert(0 && "pos_is_ok: Castling");
        }
    return true;
}

}  // namespace Stockfish
