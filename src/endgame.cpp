/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
*/

#include <algorithm>
#include <cassert>

#include "bitboard.h"
#include "endgame.h"
#include "movegen.h"

namespace Stockfish {

namespace {

inline int push_to_edge(Square s) {
    int rd = std::min(int(rank_of(s)), int(RANK_8 - rank_of(s)));
    int fd = edge_distance(file_of(s));
    return 90 - (7 * fd * fd / 2 + 7 * rd * rd / 2);
}

inline int push_to_corner(Square s) { return std::abs(7 - rank_of(s) - file_of(s)); }

inline int push_close(Square s1, Square s2) { return 140 - 20 * distance(s1, s2); }
inline int push_away(Square s1, Square s2) { return 120 - push_close(s1, s2); }

#ifndef NDEBUG
bool verify_material(const Position& pos, Color c, Value npm, int pawnsCnt) {
    return pos.non_pawn_material(c) == npm && pos.count<PAWN>(c) == pawnsCnt;
}
#endif

Square normalize(const Position& pos, Color strongSide, Square sq) {

    assert(pos.count<PAWN>(strongSide) == 1);

    if (file_of(pos.square<PAWN>(strongSide)) >= FILE_E)
        sq = flip_file(sq);

    return strongSide == WHITE ? sq : flip_rank(sq);
}

template<Color Us>
constexpr Bitboard forward_ranks_bb(Square s) {
    return Us == WHITE ? ~((Bitboard(1) << (8 * (rank_of(s) + 1))) - 1)
                       : ((Bitboard(1) << (8 * rank_of(s))) - 1);
}

inline Bitboard forward_file_bb(Color us, Square s) {
    return file_bb(s)
         & (us == WHITE ? forward_ranks_bb<WHITE>(s) : forward_ranks_bb<BLACK>(s));
}

inline Bitboard pawn_attack_span(Color us, Square s) {
    Bitboard ff = forward_file_bb(us, s);
    return (shift<EAST>(ff) & ~FileABB) | (shift<WEST>(ff) & ~FileHBB);
}

inline Bitboard passed_pawn_span(Color us, Square s) {
    return forward_file_bb(us, s) | pawn_attack_span(us, s);
}

inline Bitboard pawn_attacks_from(Color c, Square s) {
    return c == WHITE ? pawn_attacks_bb<WHITE>(square_bb(s)) : pawn_attacks_bb<BLACK>(square_bb(s));
}

inline Square frontmost_sq(Color c, Bitboard b) { return c == WHITE ? msb(b) : lsb(b); }

inline bool opposite_colors(Square s1, Square s2) { return ((int(s1) ^ int(s2)) & 1) != 0; }

inline bool pawn_passed(const Position& pos, Color c, Square s) {
    return !(passed_pawn_span(c, s) & pos.pieces(~c, PAWN));
}

}  // namespace

namespace Endgames {

std::pair<Map<Value>, Map<ScaleFactor>> maps;

void init() {

    add<KPK>("KPK");
    add<KNNK>("KNNK");
    add<KBNK>("KBNK");
    add<KRKP>("KRKP");
    add<KRKB>("KRKB");
    add<KRKN>("KRKN");
    add<KQKP>("KQKP");
    add<KQKR>("KQKR");
    add<KNNKP>("KNNKP");

    add<KRPKR>("KRPKR");
    add<KRPKB>("KRPKB");
    add<KBPKB>("KBPKB");
    add<KBPKN>("KBPKN");
    add<KBPPKB>("KBPPKB");
    add<KRPPKRP>("KRPPKRP");
}

}  // namespace Endgames

namespace Bitbases {
bool probe(Square, Square, Square, Color) { return false; }
}  // namespace Bitbases

template<>
Value Endgame<KXK>::operator()(const Position& pos) const {

    assert(verify_material(pos, weakSide, VALUE_ZERO, 0));
    assert(!pos.checkers());

    if (pos.side_to_move() == weakSide && !MoveList<LEGAL>(pos).size())
        return VALUE_DRAW;

    Square strongKing = pos.square<KING>(strongSide);
    Square weakKing   = pos.square<KING>(weakSide);

    Value result = pos.non_pawn_material(strongSide) + pos.count<PAWN>(strongSide) * PawnValueEg
                 + push_to_edge(weakKing) + push_close(strongKing, weakKing);

    if (pos.count<QUEEN>(strongSide) || pos.count<ROOK>(strongSide)
        || (pos.count<BISHOP>(strongSide) && pos.count<KNIGHT>(strongSide)))
        result = std::min(result + VALUE_KNOWN_WIN, VALUE_TB_WIN_IN_MAX_PLY - 1);

    return strongSide == pos.side_to_move() ? result : -result;
}

template<>
Value Endgame<KBNK>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, KnightValueMg + BishopValueMg, 0));
    assert(verify_material(pos, weakSide, VALUE_ZERO, 0));

    Square strongKing   = pos.square<KING>(strongSide);
    Square strongBishop = pos.square<BISHOP>(strongSide);
    Square weakKing     = pos.square<KING>(weakSide);

    Value result = (VALUE_KNOWN_WIN + 3520) + push_close(strongKing, weakKing)
                 + 420 * push_to_corner(opposite_colors(strongBishop, SQ_A1) ? flip_file(weakKing)
                                                                              : weakKing);

    return strongSide == pos.side_to_move() ? result : -result;
}

template<>
Value Endgame<KPK>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, VALUE_ZERO, 1));
    assert(verify_material(pos, weakSide, VALUE_ZERO, 0));

    Square strongKing = normalize(pos, strongSide, pos.square<KING>(strongSide));
    Square strongPawn = normalize(pos, strongSide, pos.square<PAWN>(strongSide));
    Square weakKing   = normalize(pos, strongSide, pos.square<KING>(weakSide));

    Color us = strongSide == pos.side_to_move() ? WHITE : BLACK;

    if (!Bitbases::probe(strongKing, strongPawn, weakKing, us))
        return VALUE_DRAW;

    Value result = VALUE_KNOWN_WIN + PawnValueEg + Value(rank_of(strongPawn));
    return strongSide == pos.side_to_move() ? result : -result;
}

template<>
Value Endgame<KRKP>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, RookValueMg, 0));
    assert(verify_material(pos, weakSide, VALUE_ZERO, 1));

    Square strongKing     = pos.square<KING>(strongSide);
    Square weakKing       = pos.square<KING>(weakSide);
    Square strongRook     = pos.square<ROOK>(strongSide);
    Square weakPawn       = pos.square<PAWN>(weakSide);
    Square queeningSquare = make_square(file_of(weakPawn), relative_rank(weakSide, RANK_8));
    Value  result;

    if (forward_file_bb(strongSide, strongKing) & weakPawn)
        result = RookValueEg - distance(strongKing, weakPawn);
    else if (distance(weakKing, weakPawn) >= 3 + (pos.side_to_move() == weakSide)
             && distance(weakKing, strongRook) >= 3)
        result = RookValueEg - distance(strongKing, weakPawn);
    else if (relative_rank(strongSide, weakKing) <= RANK_3 && distance(weakKing, weakPawn) == 1
             && relative_rank(strongSide, strongKing) >= RANK_4
             && distance(strongKing, weakPawn) > 2 + (pos.side_to_move() == strongSide))
        result = Value(80) - 8 * distance(strongKing, weakPawn);
    else
        result = Value(200)
               - 8 * (distance(strongKing, weakPawn + pawn_push(weakSide))
                      - distance(weakKing, weakPawn + pawn_push(weakSide))
                      - distance(weakPawn, queeningSquare));

    return strongSide == pos.side_to_move() ? result : -result;
}

template<>
Value Endgame<KRKB>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, RookValueMg, 0));
    assert(verify_material(pos, weakSide, BishopValueMg, 0));

    Value result = Value(push_to_edge(pos.square<KING>(weakSide)));
    return strongSide == pos.side_to_move() ? result : -result;
}

template<>
Value Endgame<KRKN>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, RookValueMg, 0));
    assert(verify_material(pos, weakSide, KnightValueMg, 0));

    Square weakKing   = pos.square<KING>(weakSide);
    Square weakKnight = pos.square<KNIGHT>(weakSide);
    Value  result     = Value(push_to_edge(weakKing) + push_away(weakKing, weakKnight));
    return strongSide == pos.side_to_move() ? result : -result;
}

template<>
Value Endgame<KQKP>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, QueenValueMg, 0));
    assert(verify_material(pos, weakSide, VALUE_ZERO, 1));

    Square strongKing = pos.square<KING>(strongSide);
    Square weakKing   = pos.square<KING>(weakSide);
    Square weakPawn   = pos.square<PAWN>(weakSide);

    Value result = Value(push_close(strongKing, weakKing));

    if (relative_rank(weakSide, weakPawn) != RANK_7 || distance(weakKing, weakPawn) != 1
        || ((FileBBB | FileDBB | FileEBB | FileGBB) & weakPawn))
        result += QueenValueEg - PawnValueEg;

    return strongSide == pos.side_to_move() ? result : -result;
}

template<>
Value Endgame<KQKR>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, QueenValueMg, 0));
    assert(verify_material(pos, weakSide, RookValueMg, 0));

    Square strongKing = pos.square<KING>(strongSide);
    Square weakKing   = pos.square<KING>(weakSide);

    Value result = QueenValueEg - RookValueEg + push_to_edge(weakKing) + push_close(strongKing, weakKing);

    return strongSide == pos.side_to_move() ? result : -result;
}

template<>
Value Endgame<KNNKP>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, 2 * KnightValueMg, 0));
    assert(verify_material(pos, weakSide, VALUE_ZERO, 1));

    Square weakKing = pos.square<KING>(weakSide);
    Square weakPawn = pos.square<PAWN>(weakSide);

    Value result = PawnValueEg + 2 * push_to_edge(weakKing) - 10 * relative_rank(weakSide, weakPawn);

    return strongSide == pos.side_to_move() ? result : -result;
}

template<>
Value Endgame<KNNK>::operator()(const Position&) const {
    return VALUE_DRAW;
}

template<>
ScaleFactor Endgame<KBPsK>::operator()(const Position& pos) const {

    assert(pos.non_pawn_material(strongSide) == BishopValueMg);
    assert(pos.count<PAWN>(strongSide) >= 1);

    Bitboard strongPawns = pos.pieces(strongSide, PAWN);
    Bitboard allPawns    = pos.pieces(PAWN);

    Square strongBishop = pos.square<BISHOP>(strongSide);
    Square weakKing     = pos.square<KING>(weakSide);
    Square strongKing   = pos.square<KING>(strongSide);

    if (!(strongPawns & ~FileABB) || !(strongPawns & ~FileHBB))
    {
        Square queeningSquare =
          relative_square(strongSide, make_square(file_of(lsb(strongPawns)), RANK_8));

        if (opposite_colors(queeningSquare, strongBishop) && distance(queeningSquare, weakKing) <= 1)
            return SCALE_FACTOR_DRAW;
    }

    if ((!(allPawns & ~FileBBB) || !(allPawns & ~FileGBB)) && pos.non_pawn_material(weakSide) == 0
        && pos.count<PAWN>(weakSide) >= 1)
    {
        Square weakPawn = frontmost_sq(strongSide, pos.pieces(weakSide, PAWN));

        if (relative_rank(strongSide, weakPawn) == RANK_7
            && (strongPawns & (weakPawn + pawn_push(weakSide)))
            && (opposite_colors(strongBishop, weakPawn) || !more_than_one(strongPawns)))
        {
            int strongKingDist = distance(weakPawn, strongKing);
            int weakKingDist   = distance(weakPawn, weakKing);

            if (relative_rank(strongSide, weakKing) >= RANK_7 && weakKingDist <= 2
                && weakKingDist <= strongKingDist)
                return SCALE_FACTOR_DRAW;
        }
    }

    return SCALE_FACTOR_NONE;
}

template<>
ScaleFactor Endgame<KQKRPs>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, QueenValueMg, 0));
    assert(pos.count<ROOK>(weakSide) == 1);
    assert(pos.count<PAWN>(weakSide) >= 1);

    Square strongKing = pos.square<KING>(strongSide);
    Square weakKing   = pos.square<KING>(weakSide);
    Square weakRook   = pos.square<ROOK>(weakSide);

    if (relative_rank(weakSide, weakKing) <= RANK_2 && relative_rank(weakSide, strongKing) >= RANK_4
        && relative_rank(weakSide, weakRook) == RANK_3
        && (pos.pieces(weakSide, PAWN) & attacks_bb<KING>(weakKing)
            & pawn_attacks_from(strongSide, weakRook)))
        return SCALE_FACTOR_DRAW;

    return SCALE_FACTOR_NONE;
}

template<>
ScaleFactor Endgame<KRPKR>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, RookValueMg, 1));
    assert(verify_material(pos, weakSide, RookValueMg, 0));

    Square strongKing = normalize(pos, strongSide, pos.square<KING>(strongSide));
    Square strongRook = normalize(pos, strongSide, pos.square<ROOK>(strongSide));
    Square strongPawn = normalize(pos, strongSide, pos.square<PAWN>(strongSide));
    Square weakKing   = normalize(pos, strongSide, pos.square<KING>(weakSide));
    Square weakRook   = normalize(pos, strongSide, pos.square<ROOK>(weakSide));

    File   pawnFile       = file_of(strongPawn);
    Rank   pawnRank       = rank_of(strongPawn);
    Square queeningSquare = make_square(pawnFile, RANK_8);
    int    tempo          = (pos.side_to_move() == strongSide);

    if (pawnRank <= RANK_5 && distance(weakKing, queeningSquare) <= 1 && strongKing <= SQ_H5
        && (rank_of(weakRook) == RANK_6
            || (pawnRank <= RANK_3 && rank_of(strongRook) != RANK_6)))
        return SCALE_FACTOR_DRAW;

    if (pawnRank == RANK_6 && distance(weakKing, queeningSquare) <= 1
        && rank_of(strongKing) + tempo <= RANK_6
        && (rank_of(weakRook) == RANK_1 || (!tempo && distance<File>(weakRook, strongPawn) >= 3)))
        return SCALE_FACTOR_DRAW;

    if (pawnRank >= RANK_6 && weakKing == queeningSquare && rank_of(weakRook) == RANK_1
        && (!tempo || distance(strongKing, strongPawn) >= 2))
        return SCALE_FACTOR_DRAW;

    if (strongPawn == SQ_A7 && strongRook == SQ_A8 && (weakKing == SQ_H7 || weakKing == SQ_G7)
        && file_of(weakRook) == FILE_A
        && (rank_of(weakRook) <= RANK_3 || file_of(strongKing) >= FILE_D
            || rank_of(strongKing) <= RANK_5))
        return SCALE_FACTOR_DRAW;

    if (pawnRank <= RANK_5 && weakKing == strongPawn + NORTH
        && distance(strongKing, strongPawn) - tempo >= 2
        && distance(strongKing, weakRook) - tempo >= 2)
        return SCALE_FACTOR_DRAW;

    if (pawnRank == RANK_7 && pawnFile != FILE_A && file_of(strongRook) == pawnFile
        && strongRook != queeningSquare
        && (distance(strongKing, queeningSquare) < distance(weakKing, queeningSquare) - 2 + tempo)
        && (distance(strongKing, queeningSquare) < distance(weakKing, strongRook) + tempo))
        return ScaleFactor(SCALE_FACTOR_MAX - 2 * distance(strongKing, queeningSquare));

    if (pawnFile != FILE_A && file_of(strongRook) == pawnFile && strongRook < strongPawn
        && (distance(strongKing, queeningSquare) < distance(weakKing, queeningSquare) - 2 + tempo)
        && (distance(strongKing, strongPawn + NORTH)
            < distance(weakKing, strongPawn + NORTH) - 2 + tempo)
        && (distance(weakKing, strongRook) + tempo >= 3
            || (distance(strongKing, queeningSquare) < distance(weakKing, strongRook) + tempo
                && (distance(strongKing, strongPawn + NORTH)
                    < distance(weakKing, strongPawn) + tempo))))
        return ScaleFactor(SCALE_FACTOR_MAX - 8 * distance(strongPawn, queeningSquare)
                           - 2 * distance(strongKing, queeningSquare));

    if (pawnRank <= RANK_4 && weakKing > strongPawn)
    {
        if (file_of(weakKing) == file_of(strongPawn))
            return ScaleFactor(10);
        if (distance<File>(weakKing, strongPawn) == 1 && distance(strongKing, weakKing) > 2)
            return ScaleFactor(24 - 2 * distance(strongKing, weakKing));
    }
    return SCALE_FACTOR_NONE;
}

template<>
ScaleFactor Endgame<KRPKB>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, RookValueMg, 1));
    assert(verify_material(pos, weakSide, BishopValueMg, 0));

    if (pos.pieces(PAWN) & (FileABB | FileHBB))
    {
        Square weakKing    = pos.square<KING>(weakSide);
        Square weakBishop  = pos.square<BISHOP>(weakSide);
        Square strongKing  = pos.square<KING>(strongSide);
        Square strongPawn  = pos.square<PAWN>(strongSide);
        Rank   pawnRank    = relative_rank(strongSide, strongPawn);
        Direction push     = pawn_push(strongSide);

        if (pawnRank == RANK_5 && !opposite_colors(weakBishop, strongPawn))
        {
            int d = distance(strongPawn + 3 * push, weakKing);

            if (d <= 2 && !(d == 0 && weakKing == strongKing + 2 * push))
                return ScaleFactor(24);
            else
                return ScaleFactor(48);
        }

        if (pawnRank == RANK_6 && distance(strongPawn + 2 * push, weakKing) <= 1
            && (attacks_bb<BISHOP>(weakBishop) & (strongPawn + push))
            && distance<File>(weakBishop, strongPawn) >= 2)
            return ScaleFactor(8);
    }

    return SCALE_FACTOR_NONE;
}

template<>
ScaleFactor Endgame<KRPPKRP>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, RookValueMg, 2));
    assert(verify_material(pos, weakSide, RookValueMg, 1));

    Square strongPawn1 = lsb(pos.pieces(strongSide, PAWN));
    Square strongPawn2 = msb(pos.pieces(strongSide, PAWN));
    Square weakKing    = pos.square<KING>(weakSide);

    if (pawn_passed(pos, strongSide, strongPawn1) || pawn_passed(pos, strongSide, strongPawn2))
        return SCALE_FACTOR_NONE;

    Rank pawnRank = std::max(relative_rank(strongSide, strongPawn1), relative_rank(strongSide, strongPawn2));

    if (distance<File>(weakKing, strongPawn1) <= 1 && distance<File>(weakKing, strongPawn2) <= 1
        && relative_rank(strongSide, weakKing) > pawnRank)
        return ScaleFactor(7 * pawnRank);

    return SCALE_FACTOR_NONE;
}

template<>
ScaleFactor Endgame<KPsK>::operator()(const Position& pos) const {

    assert(pos.non_pawn_material(strongSide) == VALUE_ZERO);
    assert(pos.count<PAWN>(strongSide) >= 2);
    assert(verify_material(pos, weakSide, VALUE_ZERO, 0));

    Square weakKing    = pos.square<KING>(weakSide);
    Bitboard strongPawns = pos.pieces(strongSide, PAWN);

    if (!(strongPawns & ~(FileABB | FileHBB))
        && !(strongPawns & ~passed_pawn_span(weakSide, weakKing)))
        return SCALE_FACTOR_DRAW;

    return SCALE_FACTOR_NONE;
}

template<>
ScaleFactor Endgame<KBPKB>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, BishopValueMg, 1));
    assert(verify_material(pos, weakSide, BishopValueMg, 0));

    Square strongPawn   = pos.square<PAWN>(strongSide);
    Square strongBishop = pos.square<BISHOP>(strongSide);
    Square weakBishop   = pos.square<BISHOP>(weakSide);
    Square weakKing     = pos.square<KING>(weakSide);

    if ((forward_file_bb(strongSide, strongPawn) & weakKing)
        && (opposite_colors(weakKing, strongBishop)
            || relative_rank(strongSide, weakKing) <= RANK_6))
        return SCALE_FACTOR_DRAW;

    if (opposite_colors(strongBishop, weakBishop))
        return SCALE_FACTOR_DRAW;

    return SCALE_FACTOR_NONE;
}

template<>
ScaleFactor Endgame<KBPPKB>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, BishopValueMg, 2));
    assert(verify_material(pos, weakSide, BishopValueMg, 0));

    Square strongBishop = pos.square<BISHOP>(strongSide);
    Square weakBishop   = pos.square<BISHOP>(weakSide);

    if (!opposite_colors(strongBishop, weakBishop))
        return SCALE_FACTOR_NONE;

    Square weakKing    = pos.square<KING>(weakSide);
    Square strongPawn1 = lsb(pos.pieces(strongSide, PAWN));
    Square strongPawn2 = msb(pos.pieces(strongSide, PAWN));
    Square blockSq1, blockSq2;

    if (relative_rank(strongSide, strongPawn1) > relative_rank(strongSide, strongPawn2))
    {
        blockSq1 = strongPawn1 + pawn_push(strongSide);
        blockSq2 = make_square(file_of(strongPawn2), rank_of(strongPawn1));
    }
    else
    {
        blockSq1 = strongPawn2 + pawn_push(strongSide);
        blockSq2 = make_square(file_of(strongPawn1), rank_of(strongPawn2));
    }

    switch (distance<File>(strongPawn1, strongPawn2))
    {
    case 0:
        if (file_of(weakKing) == file_of(blockSq1)
            && relative_rank(strongSide, weakKing) >= relative_rank(strongSide, blockSq1)
            && opposite_colors(weakKing, strongBishop))
            return SCALE_FACTOR_DRAW;
        else
            return SCALE_FACTOR_NONE;

    case 1:
        if (weakKing == blockSq1 && opposite_colors(weakKing, strongBishop)
            && (weakBishop == blockSq2
                || (attacks_bb<BISHOP>(blockSq2, pos.pieces()) & pos.pieces(weakSide, BISHOP))
                || distance<Rank>(strongPawn1, strongPawn2) >= 2))
            return SCALE_FACTOR_DRAW;

        else if (weakKing == blockSq2 && opposite_colors(weakKing, strongBishop)
                 && (weakBishop == blockSq1
                     || (attacks_bb<BISHOP>(blockSq1, pos.pieces())
                         & pos.pieces(weakSide, BISHOP))))
            return SCALE_FACTOR_DRAW;
        else
            return SCALE_FACTOR_NONE;

    default: return SCALE_FACTOR_NONE;
    }
}

template<>
ScaleFactor Endgame<KBPKN>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, BishopValueMg, 1));
    assert(verify_material(pos, weakSide, KnightValueMg, 0));

    Square strongPawn   = pos.square<PAWN>(strongSide);
    Square strongBishop = pos.square<BISHOP>(strongSide);
    Square weakKing     = pos.square<KING>(weakSide);

    if (file_of(weakKing) == file_of(strongPawn)
        && relative_rank(strongSide, strongPawn) < relative_rank(strongSide, weakKing)
        && (opposite_colors(weakKing, strongBishop)
            || relative_rank(strongSide, weakKing) <= RANK_6))
        return SCALE_FACTOR_DRAW;

    return SCALE_FACTOR_NONE;
}

template<>
ScaleFactor Endgame<KPKP>::operator()(const Position& pos) const {

    assert(verify_material(pos, strongSide, VALUE_ZERO, 1));
    assert(verify_material(pos, weakSide, VALUE_ZERO, 1));

    Square strongKing = normalize(pos, strongSide, pos.square<KING>(strongSide));
    Square weakKing   = normalize(pos, strongSide, pos.square<KING>(weakSide));
    Square strongPawn = normalize(pos, strongSide, pos.square<PAWN>(strongSide));

    Color us = strongSide == pos.side_to_move() ? WHITE : BLACK;

    if (rank_of(strongPawn) >= RANK_5 && file_of(strongPawn) != FILE_A)
        return SCALE_FACTOR_NONE;

    return Bitbases::probe(strongKing, strongPawn, weakKing, us) ? SCALE_FACTOR_NONE : SCALE_FACTOR_DRAW;
}

}  // namespace Stockfish
