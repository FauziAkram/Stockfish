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

#include "evaluate.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>

#include "bitboard.h"
#include "material.h"
#include "position.h"
#include "uci.h"
#include "pawns.h"
#include "thread.h"

using namespace std;

namespace Stockfish {

namespace Eval {

}

namespace Trace {

  enum Tracing { NO_TRACE, TRACE };

  enum Term { // The first 8 entries are reserved for PieceType
    MATERIAL = 8, IMBALANCE, MOBILITY, THREAT, PASSED, SPACE, WINNABLE, TOTAL, TERM_NB
  };

  Score scores[TERM_NB][COLOR_NB];

  static double to_cp(Value v) { return double(v) / 100; }

  static void add(int idx, Color c, Score s) {
    scores[idx][c] = s;
  }

  static void add(int idx, Score w, Score b = SCORE_ZERO) {
    scores[idx][WHITE] = w;
    scores[idx][BLACK] = b;
  }

  static std::ostream& operator<<(std::ostream& os, Score s) {
    os << std::setw(5) << to_cp(mg_value(s)) << " "
       << std::setw(5) << to_cp(eg_value(s));
    return os;
  }

  static std::ostream& operator<<(std::ostream& os, Term t) {

    if (t == MATERIAL || t == IMBALANCE || t == WINNABLE || t == TOTAL)
        os << " ----  ----"    << " | " << " ----  ----";
    else
        os << scores[t][WHITE] << " | " << scores[t][BLACK];

    os << " | " << scores[t][WHITE] - scores[t][BLACK] << " |\n";
    return os;
  }
}

using namespace Trace;

namespace {

  // Threshold for lazy and space evaluation
  constexpr Value LazyThreshold1    =  Value(3622);
  constexpr Value LazyThreshold2    =  Value(1962);
  constexpr Value SpaceThreshold    =  Value(11551);

  // KingAttackWeights[PieceType] contains king attack weights by piece type
  constexpr int KingAttackWeights[PIECE_TYPE_NB] = { 0, 0, 76, 46, 45, 14 };

  // SafeCheck[PieceType][single/multiple] contains safe check bonus by piece type,
  // higher if multiple safe checks are possible for that piece type.
  constexpr int SafeCheck[][2] = {
      {}, {}, {805, 1292}, {650, 984}, {1071, 1886}, {730, 1128}
  };

#define S(mg, eg) make_score(mg, eg)

  // MobilityBonus[PieceType-2][attacked] contains bonuses for middle and end game,
  // indexed by piece type and number of attacked squares in the mobility area.
  constexpr Score MobilityBonus[][32] = {
    { S(-62,-79), S(-53,-57), S(-12,-31), S( -3,-17), S(  3,  7), S( 12, 13), // Knight
      S( 21, 16), S( 28, 21), S( 37, 26) },
    { S(-47,-59), S(-20,-25), S( 14, -8), S( 29, 12), S( 39, 21), S( 53, 40), // Bishop
      S( 53, 56), S( 60, 58), S( 62, 65), S( 69, 72), S( 78, 78), S( 83, 87),
      S( 91, 88), S( 96, 98) },
    { S(-60,-82), S(-24,-15), S(  0, 17) ,S(  3, 43), S(  4, 72), S( 14,100), // Rook
      S( 20,102), S( 30,122), S( 41,133), S(41 ,139), S( 41,153), S( 45,160),
      S( 57,165), S( 58,170), S( 67,175) },
    { S(-29,-49), S(-16,-29), S( -8, -8), S( -8, 17), S( 18, 39), S( 25, 54), // Queen
      S( 23, 59), S( 37, 73), S( 41, 76), S( 54, 95), S( 65, 95) ,S( 68,101),
      S( 69,124), S( 70,128), S( 70,132), S( 70,133) ,S( 71,136), S( 72,140),
      S( 74,147), S( 76,149), S( 90,153), S(104,169), S(105,171), S(106,171),
      S(112,178), S(114,185), S(114,187), S(119,221) }
  };

  // BishopPawns[distance from edge] contains a file-dependent penalty for pawns on
  // squares of the same color as our bishop.
  constexpr Score BishopPawns[int(FILE_NB) / 2] = {
    S(3, 8), S(3, 9), S(2, 7), S(3, 7)
  };

  // KingProtector[knight/bishop] contains penalty for each distance unit to own king
  constexpr Score KingProtector[] = { S(9, 9), S(7, 9) };

  // Outpost[knight/bishop] contains bonuses for each knight or bishop occupying a
  // pawn protected square on rank 4 to 6 which is also safe from a pawn attack.
  constexpr Score Outpost[] = { S(54, 34), S(31, 25) };

  // PassedRank[Rank] contains a bonus according to the rank of a passed pawn
  constexpr Score PassedRank[RANK_NB] = {
    S(0, 0), S(2, 38), S(15, 36), S(22, 50), S(64, 81), S(166, 184), S(284, 269)
  };

  constexpr Score RookOnClosedFile = S(10, 5);
  constexpr Score RookOnOpenFile[] = { S(18, 8), S(49, 26) };

  // ThreatByMinor/ByRook[attacked PieceType] contains bonuses according to
  // which piece type attacks which one. Attacks on lesser pieces which are
  // pawn-defended are not considered.
  constexpr Score ThreatByMinor[PIECE_TYPE_NB] = {
    S(0, 0), S(6, 37), S(64, 50), S(82, 57), S(103, 130), S(81, 163)
  };

  constexpr Score ThreatByRook[PIECE_TYPE_NB] = {
    S(0, 0), S(3, 44), S(36, 71), S(44, 59), S(0, 39), S(60, 39)
  };

  constexpr Value CorneredBishop = Value(50);

  // Assorted bonuses and penalties
  constexpr Score UncontestedOutpost  = S(  0, 10);
  constexpr Score BishopOnKingRing    = S( 24,  0);
  constexpr Score BishopXRayPawns     = S(  4,  5);
  constexpr Score FlankAttacks        = S(  8,  0);
  constexpr Score Hanging             = S( 72, 40);
  constexpr Score KnightOnQueen       = S( 16, 11);
  constexpr Score LongDiagonalBishop  = S( 45,  0);
  constexpr Score MinorBehindPawn     = S( 18,  3);
  constexpr Score PassedFile          = S( 13,  8);
  constexpr Score PawnlessFlank       = S( 19, 97);
  constexpr Score ReachableOutpost    = S( 33, 19);
  constexpr Score RestrictedPiece     = S(  6,  7);
  constexpr Score RookOnKingRing      = S( 16,  0);
  constexpr Score SliderOnQueen       = S( 62, 21);
  constexpr Score ThreatByKing        = S( 24, 87);
  constexpr Score ThreatByPawnPush    = S( 48, 39);
  constexpr Score ThreatBySafePawn    = S(167, 99);
  constexpr Score TrappedRook         = S( 55, 13);
  constexpr Score WeakQueenProtection = S( 14,  0);
  constexpr Score WeakQueen           = S( 57, 19);


#undef S

  // Evaluation class computes and stores attacks tables and other working data
  template<Tracing T>
  class Evaluation {

  public:
    Evaluation() = delete;
    explicit Evaluation(const Position& p) : pos(p) {}
    Evaluation& operator=(const Evaluation&) = delete;
    Value value();

  private:
    template<Color Us> void initialize();
    template<Color Us, PieceType Pt> Score pieces();
    template<Color Us> Score king() const;
    template<Color Us> Score threats() const;
    template<Color Us> Score passed() const;
    template<Color Us> Score space() const;
    Value winnable(Score score) const;

    const Position& pos;
    Material::Entry* me;
    Pawns::Entry* pe;
    Bitboard mobilityArea[COLOR_NB];
    Score mobility[COLOR_NB] = { SCORE_ZERO, SCORE_ZERO };

    // attackedBy[color][piece type] is a bitboard representing all squares
    // attacked by a given color and piece type. Special "piece types" which
    // is also calculated is ALL_PIECES.
    Bitboard attackedBy[COLOR_NB][PIECE_TYPE_NB];

    // attackedBy2[color] are the squares attacked by at least 2 units of a given
    // color, including x-rays. But diagonal x-rays through pawns are not computed.
    Bitboard attackedBy2[COLOR_NB];

    // kingRing[color] are the squares adjacent to the king plus some other
    // very near squares, depending on king position.
    Bitboard kingRing[COLOR_NB];

    // kingAttackersCount[color] is the number of pieces of the given color
    // which attack a square in the kingRing of the enemy king.
    int kingAttackersCount[COLOR_NB];

    // kingAttackersWeight[color] is the sum of the "weights" of the pieces of
    // the given color which attack a square in the kingRing of the enemy king.
    // The weights of the individual piece types are given by the elements in
    // the KingAttackWeights array.
    int kingAttackersWeight[COLOR_NB];

    // kingAttacksCount[color] is the number of attacks by the given color to
    // squares directly adjacent to the enemy king. Pieces which attack more
    // than one square are counted multiple times. For instance, if there is
    // a white knight on g5 and black's king is on g8, this white knight adds 2
    // to kingAttacksCount[WHITE].
    int kingAttacksCount[COLOR_NB];
  };


  // Evaluation::initialize() computes king and pawn attacks, and the king ring
  // bitboard for a given color. This is done at the beginning of the evaluation.

  template<Tracing T> template<Color Us>
  void Evaluation<T>::initialize() {

    constexpr Color     Them = ~Us;
    constexpr Direction Up   = pawn_push(Us);
    constexpr Direction Down = -Up;
    constexpr Bitboard LowRanks = (Us == WHITE ? Rank2BB | Rank3BB : Rank7BB | Rank6BB);

    const Square ksq = pos.square<KING>(Us);

    Bitboard dblAttackByPawn = pawn_double_attacks_bb<Us>(pos.pieces(Us, PAWN));

    // Find our pawns that are blocked or on the first two ranks
    Bitboard b = pos.pieces(Us, PAWN) & (shift<Down>(pos.pieces()) | LowRanks);

    // Squares occupied by those pawns, by our king or queen, by blockers to attacks on our king
    // or controlled by enemy pawns are excluded from the mobility area.
    mobilityArea[Us] = ~(b | pos.pieces(Us, KING, QUEEN) | pos.blockers_for_king(Us) | pe->pawn_attacks(Them));

    // Initialize attackedBy[] for king and pawns
    attackedBy[Us][KING] = attacks_bb<KING>(ksq);
    attackedBy[Us][PAWN] = pe->pawn_attacks(Us);
    attackedBy[Us][ALL_PIECES] = attackedBy[Us][KING] | attackedBy[Us][PAWN];
    attackedBy2[Us] = dblAttackByPawn | (attackedBy[Us][KING] & attackedBy[Us][PAWN]);

    // Init our king safety tables
    Square s = make_square(std::clamp(file_of(ksq), FILE_B, FILE_G),
                           std::clamp(rank_of(ksq), RANK_2, RANK_7));
    kingRing[Us] = attacks_bb<KING>(s) | s;

    kingAttackersCount[Them] = popcount(kingRing[Us] & pe->pawn_attacks(Them));
    kingAttacksCount[Them] = kingAttackersWeight[Them] = 0;

    // Remove from kingRing[] the squares defended by two pawns
    kingRing[Us] &= ~dblAttackByPawn;
  }


  // Evaluation::pieces() scores pieces of a given color and type

  template<Tracing T> template<Color Us, PieceType Pt>
  Score Evaluation<T>::pieces() {

    constexpr Color Them = ~Us;
    [[maybe_unused]] constexpr Direction Down = -pawn_push(Us);
    [[maybe_unused]] constexpr Bitboard OutpostRanks = (Us == WHITE ? Rank4BB | Rank5BB | Rank6BB
                                                                    : Rank5BB | Rank4BB | Rank3BB);
    Bitboard b1 = pos.pieces(Us, Pt);
    Bitboard b, bb;
    Score score = SCORE_ZERO;

    attackedBy[Us][Pt] = 0;

    while (b1)
    {
        Square s = pop_lsb(b1);

        // Find attacked squares, including x-ray attacks for bishops and rooks
        b = Pt == BISHOP ? attacks_bb<BISHOP>(s, pos.pieces() ^ pos.pieces(QUEEN))
          : Pt ==   ROOK ? attacks_bb<  ROOK>(s, pos.pieces() ^ pos.pieces(QUEEN) ^ pos.pieces(Us, ROOK))
                         : attacks_bb<Pt>(s, pos.pieces());

        if (pos.blockers_for_king(Us) & s)
            b &= line_bb(pos.square<KING>(Us), s);

        attackedBy2[Us] |= attackedBy[Us][ALL_PIECES] & b;
        attackedBy[Us][Pt] |= b;
        attackedBy[Us][ALL_PIECES] |= b;

        if (b & kingRing[Them])
        {
            kingAttackersCount[Us]++;
            kingAttackersWeight[Us] += KingAttackWeights[Pt];
            kingAttacksCount[Us] += popcount(b & attacks_bb<KING>(pos.square<KING>(Them)));
        }
        else if (Pt == ROOK && (file_bb(s) & kingRing[Them]))
            score += RookOnKingRing;
        else if (Pt == BISHOP && (attacks_bb<BISHOP>(s, pos.pieces(PAWN)) & kingRing[Them]))
            score += BishopOnKingRing;

        int mob = popcount(b & mobilityArea[Us]);
        mobility[Us] += MobilityBonus[Pt - 2][mob];

        if constexpr (Pt == BISHOP || Pt == KNIGHT)
        {
            // Bonus if the piece is on an outpost square or can reach one
            // Bonus for knights (UncontestedOutpost) if few relevant targets
            bb = OutpostRanks & (attackedBy[Us][PAWN] | shift<Down>(pos.pieces(PAWN)))
                              & ~pe->pawn_attacks_span(Them);
            Bitboard targets = pos.pieces(Them) & ~pos.pieces(PAWN);

            if (   Pt == KNIGHT
                && (bb & s & ~CenterFiles) // on a side outpost
                && !(b & targets)        // no relevant attacks
                && (!more_than_one(targets & (s & QueenSide ? QueenSide : KingSide))))
                score += UncontestedOutpost * popcount(pos.pieces(PAWN) & (s & QueenSide ? QueenSide : KingSide));
            else if (bb & s)
                score += Outpost[Pt == BISHOP];
            else if (Pt == KNIGHT && (bb & b & ~pos.pieces(Us)))
                score += ReachableOutpost;

            // Bonus for a knight or bishop shielded by pawn
            if (shift<Down>(pos.pieces(PAWN)) & s)
                score += MinorBehindPawn;

            // Penalty if the piece is far from the king
            score -= KingProtector[Pt == BISHOP] * distance(pos.square<KING>(Us), s);

            if constexpr (Pt == BISHOP)
            {
                // Penalty according to the number of our pawns on the same color square as the
                // bishop, bigger when the center files are blocked with pawns and smaller
                // when the bishop is outside the pawn chain.
                Bitboard blocked = pos.pieces(Us, PAWN) & shift<Down>(pos.pieces());

                score -= BishopPawns[edge_distance(file_of(s))] * pos.pawns_on_same_color_squares(Us, s)
                                     * (!(attackedBy[Us][PAWN] & s) + popcount(blocked & CenterFiles));

                // Penalty for all enemy pawns x-rayed
                score -= BishopXRayPawns * popcount(attacks_bb<BISHOP>(s) & pos.pieces(Them, PAWN));

                // Bonus for bishop on a long diagonal which can "see" both center squares
                if (more_than_one(attacks_bb<BISHOP>(s, pos.pieces(PAWN)) & Center))
                    score += LongDiagonalBishop;

                // An important Chess960 pattern: a cornered bishop blocked by a friendly
                // pawn diagonally in front of it is a very serious problem, especially
                // when that pawn is also blocked.
                if (   pos.is_chess960()
                    && (s == relative_square(Us, SQ_A1) || s == relative_square(Us, SQ_H1)))
                {
                    Direction d = pawn_push(Us) + (file_of(s) == FILE_A ? EAST : WEST);
                    if (pos.piece_on(s + d) == make_piece(Us, PAWN))
                        score -= !pos.empty(s + d + pawn_push(Us)) ? 4 * make_score(CorneredBishop, CorneredBishop)
                                                                   : 3 * make_score(CorneredBishop, CorneredBishop);
                }
            }
        }

        if constexpr (Pt == ROOK)
        {
            // Bonuses for rook on a (semi-)open or closed file
            if (pos.is_on_semiopen_file(Us, s))
            {
                score += RookOnOpenFile[pos.is_on_semiopen_file(Them, s)];
            }
            else
            {
                // If our pawn on this file is blocked, increase penalty
                if ( pos.pieces(Us, PAWN)
                   & shift<Down>(pos.pieces())
                   & file_bb(s))
                {
                    score -= RookOnClosedFile;
                }

                // Penalty when trapped by the king, even more if the king cannot castle
                if (mob <= 3)
                {
                    File kf = file_of(pos.square<KING>(Us));
                    if ((kf < FILE_E) == (file_of(s) < kf))
                        score -= TrappedRook * (1 + !pos.can_castle(Us & ANY_CASTLING));
                }
            }
        }

        if constexpr (Pt == QUEEN)
        {
            // Penalty if any relative pin or discovered attack against the queen
            Bitboard queenPinners;
            if (pos.slider_blockers(pos.pieces(Them, ROOK, BISHOP), s, queenPinners))
                score -= WeakQueen;
        }
    }
    if constexpr (T)
        Trace::add(Pt, Us, score);

    return score;
  }

  template<Tracing T>
  Value Evaluation<T>::value() {

    assert(!pos.checkers());

    // Probe the material hash table
    me = Material::probe(pos);

    // If we have a specialized evaluation function for the current material
    // configuration, call it and return.
    if (me->specialized_eval_exists())
        return me->evaluate(pos);

    // Initialize score by reading the incrementally updated scores included in
    // the position object (material + piece square tables) and the material
    // imbalance. Score is computed internally from the white point of view.
    Score score = pos.psq_score() + me->imbalance();

    // Probe the pawn hash table
    pe = Pawns::probe(pos);
    score += pe->pawn_score(WHITE) - pe->pawn_score(BLACK);

    // Early exit if score is high
    auto lazy_skip = [&](Value lazyThreshold) {
        return abs(mg_value(score) + eg_value(score)) >   lazyThreshold
                                                        + std::abs(pos.this_thread()->bestValue) * 5 / 4
                                                        + pos.non_pawn_material() / 32;
    };

    if (lazy_skip(LazyThreshold1))
        goto make_v;

    // Main evaluation begins here
    initialize<WHITE>();
    initialize<BLACK>();

    // Pieces evaluated first (also populates attackedBy, attackedBy2).
    // Note that the order of evaluation of the terms is left unspecified.
    score +=  pieces<WHITE, KNIGHT>() - pieces<BLACK, KNIGHT>()
            + pieces<WHITE, BISHOP>() - pieces<BLACK, BISHOP>()
            + pieces<WHITE, ROOK  >() - pieces<BLACK, ROOK  >()
            + pieces<WHITE, QUEEN >() - pieces<BLACK, QUEEN >();

    score += mobility[WHITE] - mobility[BLACK];

    // More complex interactions that require fully populated attack bitboards
    score +=  king<   WHITE>() - king<   BLACK>()
            + passed< WHITE>() - passed< BLACK>();

    if (lazy_skip(LazyThreshold2))
        goto make_v;

    score +=  threats<WHITE>() - threats<BLACK>()
            + space<  WHITE>() - space<  BLACK>();

make_v:
    // Derive single value from mg and eg parts of score
    Value v = winnable(score);

    // In case of tracing add all remaining individual evaluation terms
    if constexpr (T)
    {
        Trace::add(MATERIAL, pos.psq_score());
        Trace::add(IMBALANCE, me->imbalance());
        Trace::add(PAWN, pe->pawn_score(WHITE), pe->pawn_score(BLACK));
        Trace::add(MOBILITY, mobility[WHITE], mobility[BLACK]);
    }

    // Evaluation grain
    v = (v / 16) * 16;

    // Side to move point of view
    v = (pos.side_to_move() == WHITE ? v : -v);

    return v;
  }

} // namespace

Value Eval::evaluate(const Position& pos) {

  assert(!pos.checkers());

  Value v = Evaluation<NO_TRACE>(pos).value();

  // Damp down the evaluation linearly when shuffling
  v = v * (200 - pos.rule50_count()) / 214;

  // Guarantee evaluation does not hit the tablebase range
  v = std::clamp(v, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);

  return v;
}

std::string Eval::trace(Position& pos) {

  if (pos.checkers())
      return "Final evaluation: none (in check)";

  std::stringstream ss;
  ss << std::showpoint << std::noshowpos << std::fixed << std::setprecision(2);

  Value v;

  std::memset(scores, 0, sizeof(scores));

  // Reset any global variable used in eval
  pos.this_thread()->bestValue       = VALUE_ZERO;

  v = Evaluation<TRACE>(pos).value();

  ss << std::showpoint << std::noshowpos << std::fixed << std::setprecision(2)
     << " Contributing terms for the classical eval:\n"
     << "+------------+-------------+-------------+-------------+\n"
     << "|    Term    |    White    |    Black    |    Total    |\n"
     << "|            |   MG    EG  |   MG    EG  |   MG    EG  |\n"
     << "+------------+-------------+-------------+-------------+\n"
     << "|   Material | " << Term(MATERIAL)
     << "|  Imbalance | " << Term(IMBALANCE)
     << "|      Pawns | " << Term(PAWN)
     << "|    Knights | " << Term(KNIGHT)
     << "|    Bishops | " << Term(BISHOP)
     << "|      Rooks | " << Term(ROOK)
     << "|     Queens | " << Term(QUEEN)
     << "|   Mobility | " << Term(MOBILITY)
     << "|King safety | " << Term(KING)
     << "|    Threats | " << Term(THREAT)
     << "|     Passed | " << Term(PASSED)
     << "|      Space | " << Term(SPACE)
     << "|   Winnable | " << Term(WINNABLE)
     << "+------------+-------------+-------------+-------------+\n"
     << "|      Total | " << Term(TOTAL)
     << "+------------+-------------+-------------+-------------+\n";

  ss << std::showpoint << std::showpos << std::fixed << std::setprecision(2) << std::setw(15);

  v = pos.side_to_move() == WHITE ? v : -v;
  ss << "\nFinal evaluation       " << to_cp(v) << " (white side)\n";

  return ss.str();
}

}  // namespace Stockfish
