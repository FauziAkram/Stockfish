/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
*/

#include "pawns.h"

#include <algorithm>
#include <cassert>

#include "bitboard.h"
#include "position.h"

namespace Stockfish {

namespace {

#define V Value
#define S(mg, eg) make_score(mg, eg)

constexpr EvalScore Backward      = S(6, 19);
constexpr EvalScore Doubled       = S(11, 51);
constexpr EvalScore DoubledEarly  = S(17, 7);
constexpr EvalScore Isolated      = S(1, 20);
constexpr EvalScore WeakLever     = S(2, 57);
constexpr EvalScore WeakUnopposed = S(15, 18);

constexpr EvalScore BlockedPawn[2] = {S(-19, -8), S(-7, 3)};

constexpr EvalScore BlockedStorm[RANK_NB] = {
  S(0, 0), S(0, 0), S(64, 75), S(-3, 14), S(-12, 19), S(-7, 4), S(-10, 5)};

constexpr int Connected[RANK_NB] = {0, 3, 7, 7, 15, 54, 86};

constexpr Value ShelterStrength[int(FILE_NB) / 2][RANK_NB] = {
  {V(-2), V(85), V(95), V(53), V(39), V(23), V(25)},
  {V(-55), V(64), V(32), V(-55), V(-30), V(-11), V(-61)},
  {V(-11), V(75), V(19), V(-6), V(26), V(9), V(-47)},
  {V(-41), V(-11), V(-27), V(-58), V(-42), V(-66), V(-163)}};

constexpr Value UnblockedStorm[int(FILE_NB) / 2][RANK_NB] = {
  {V(94), V(-280), V(-170), V(90), V(59), V(47), V(53)},
  {V(43), V(-17), V(128), V(39), V(26), V(-17), V(15)},
  {V(-9), V(62), V(170), V(34), V(-5), V(-20), V(-11)},
  {V(-27), V(-19), V(106), V(10), V(2), V(-13), V(-24)}};

constexpr EvalScore KingOnFile[2][2] = {{S(-18, 11), S(-6, -3)}, {S(0, 0), S(5, -4)}};

#undef S
#undef V

template<Color Us>
constexpr Bitboard forward_ranks_bb(Square s) {
    return Us == WHITE ? ~((Bitboard(1) << (8 * (rank_of(s) + 1))) - 1)
                       : ((Bitboard(1) << (8 * rank_of(s))) - 1);
}

template<Color Us>
constexpr Bitboard forward_file_bb(Square s) {
    return file_bb(s) & forward_ranks_bb<Us>(s);
}

template<Color Us>
constexpr Bitboard pawn_attack_span(Square s) {
    Bitboard ff = forward_file_bb<Us>(s);
    return (shift<EAST>(ff) & ~FileABB) | (shift<WEST>(ff) & ~FileHBB);
}

template<Color Us>
constexpr Bitboard passed_pawn_span(Square s) {
    return forward_file_bb<Us>(s) | pawn_attack_span<Us>(s);
}

constexpr Bitboard adjacent_files_bb(Square s) {
    return ((file_bb(s) & ~FileABB) >> 1) | ((file_bb(s) & ~FileHBB) << 1);
}

template<Color C>
constexpr Bitboard pawn_double_attacks_bb(Bitboard b) {
    return C == WHITE ? shift<NORTH_WEST>(b) & shift<NORTH_EAST>(b)
                      : shift<SOUTH_WEST>(b) & shift<SOUTH_EAST>(b);
}

template<Color Us>
Square frontmost_sq(Bitboard b) {
    return Us == WHITE ? msb(b) : lsb(b);
}

template<Color Us>
EvalScore evaluate(const Position& pos, Pawns::Entry* e) {

    constexpr Color     Them = ~Us;
    constexpr Direction Up   = pawn_push(Us);
    constexpr Direction Down = Direction(-Up);

    Bitboard neighbours, stoppers, support, phalanx, opposed;
    Bitboard lever, leverPush, blocked;
    Square   s;
    bool     backward, passed, doubled;
    EvalScore score = 0;
    Bitboard  b     = pos.pieces(Us, PAWN);

    Bitboard ourPawns   = pos.pieces(Us, PAWN);
    Bitboard theirPawns = pos.pieces(Them, PAWN);

    Bitboard doubleAttackThem = pawn_double_attacks_bb<Them>(theirPawns);

    e->passedPawns[Us]  = 0;
    e->kingSquares[Us]  = SQ_NONE;
    e->pawnAttacks[Us]  = e->pawnAttacksSpan[Us] = pawn_attacks_bb<Us>(ourPawns);
    e->blockedCount += popcount(shift<Up>(ourPawns) & (theirPawns | doubleAttackThem));

    while (b)
    {
        s = pop_lsb(b);

        assert(pos.piece_on(s) == make_piece(Us, PAWN));

        Rank r = relative_rank(Us, s);

        opposed    = theirPawns & forward_file_bb<Us>(s);
        blocked    = theirPawns & (s + Up);
        stoppers   = theirPawns & passed_pawn_span<Us>(s);
        lever      = theirPawns & pawn_attacks_bb<Us>(s);
        leverPush  = theirPawns & pawn_attacks_bb<Us>(s + Up);
        doubled    = ourPawns & (s - Up);
        neighbours = ourPawns & adjacent_files_bb(s);
        phalanx    = neighbours & rank_bb(s);
        support    = neighbours & rank_bb(s - Up);

        if (doubled)
        {
            if (!(ourPawns & shift<Down>(theirPawns | pawn_attacks_bb<Them>(theirPawns))))
                score -= DoubledEarly;
        }

        backward = !(neighbours & forward_ranks_bb<Them>(s + Up)) && (leverPush | blocked);

        if (!backward && !blocked)
            e->pawnAttacksSpan[Us] |= pawn_attack_span<Us>(s);

        passed = !(stoppers ^ lever)
              || (!(stoppers ^ leverPush) && popcount(phalanx) >= popcount(leverPush))
              || (stoppers == blocked && r >= RANK_5
                  && (shift<Up>(support) & ~(theirPawns | doubleAttackThem)));

        passed &= !(forward_file_bb<Us>(s) & ourPawns);

        if (passed)
            e->passedPawns[Us] |= s;

        if (support | phalanx)
        {
            int v = Connected[r] * (2 + bool(phalanx) - bool(opposed)) + 22 * popcount(support);
            score += make_score(v, v * (r - 2) / 4);
        }
        else if (!neighbours)
        {
            if (opposed && (ourPawns & forward_file_bb<Them>(s))
                && !(theirPawns & adjacent_files_bb(s)))
                score -= Doubled;
            else
                score -= Isolated + WeakUnopposed * !opposed;
        }
        else if (backward)
            score -= Backward + WeakUnopposed * !opposed * bool(~(FileABB | FileHBB) & s);

        if (!support)
            score -= Doubled * doubled + WeakLever * more_than_one(lever);

        if (blocked && r >= RANK_5)
            score += BlockedPawn[r - RANK_5];
    }

    return score;
}

}  // namespace

namespace Pawns {

Entry* probe(const Position& pos) {

    Key    key = pos.pawn_key();
    thread_local Table localPawnsTable;
    Entry* e   = localPawnsTable[key];

    if (e->key == key)
        return e;

    *e               = Entry{};
    e->key           = key;
    e->blockedCount  = 0;
    e->scores[WHITE] = evaluate<WHITE>(pos, e);
    e->scores[BLACK] = evaluate<BLACK>(pos, e);

    return e;
}

template<Color Us>
EvalScore Entry::evaluate_shelter(const Position& pos, Square ksq) const {

    constexpr Color Them = ~Us;

    Bitboard b          = pos.pieces(PAWN) & ~forward_ranks_bb<Them>(ksq);
    Bitboard ourPawns   = b & pos.pieces(Us) & ~pawnAttacks[Them];
    Bitboard theirPawns = b & pos.pieces(Them);

    EvalScore bonus = make_score(5, 5);

    File center = std::clamp(file_of(ksq), FILE_B, FILE_G);
    for (File f = File(center - 1); f <= File(center + 1); ++f)
    {
        b           = ourPawns & file_bb(f);
        int ourRank = b ? relative_rank(Us, frontmost_sq<Them>(b)) : 0;

        b             = theirPawns & file_bb(f);
        int theirRank = b ? relative_rank(Us, frontmost_sq<Them>(b)) : 0;

        int d = edge_distance(f);
        bonus += make_score(ShelterStrength[d][ourRank], 0);

        if (ourRank && (ourRank == theirRank - 1))
            bonus -= BlockedStorm[theirRank];
        else
            bonus -= make_score(UnblockedStorm[d][theirRank], 0);
    }

    const bool usSemiOpen   = !(pos.pieces(Us, PAWN) & file_bb(ksq));
    const bool themSemiOpen = !(pos.pieces(Them, PAWN) & file_bb(ksq));
    bonus -= KingOnFile[usSemiOpen][themSemiOpen];

    return bonus;
}

template<Color Us>
EvalScore Entry::do_king_safety(const Position& pos) {

    Square ksq = pos.square<KING>(Us);
    kingSquares[Us] = ksq;
    castlingRights[Us] = int(pos.state()->castlingRights & (Us == WHITE ? WHITE_CASTLING : BLACK_CASTLING));

    EvalScore shelter = evaluate_shelter<Us>(pos, ksq);

    if (pos.can_castle(Us & KING_SIDE))
        shelter = std::max(shelter, evaluate_shelter<Us>(pos, relative_square(Us, SQ_G1)));

    if (pos.can_castle(Us & QUEEN_SIDE))
        shelter = std::max(shelter, evaluate_shelter<Us>(pos, relative_square(Us, SQ_C1)));

    Bitboard pawns = pos.pieces(Us, PAWN);
    int      minPawnDist = 6;

    if (pawns & attacks_bb<KING>(ksq))
        minPawnDist = 1;
    else
        while (pawns)
            minPawnDist = std::min(minPawnDist, distance(ksq, pop_lsb(pawns)));

    return shelter - make_score(0, 16 * minPawnDist);
}

template EvalScore Entry::do_king_safety<WHITE>(const Position& pos);
template EvalScore Entry::do_king_safety<BLACK>(const Position& pos);

template EvalScore Entry::evaluate_shelter<WHITE>(const Position& pos, Square ksq) const;
template EvalScore Entry::evaluate_shelter<BLACK>(const Position& pos, Square ksq) const;

}  // namespace Pawns

}  // namespace Stockfish
