#include "kpk_bitbase.h"

#include <vector>
#include <bitset>
#include <cassert>

#include "bitboard.h"
#include "misc.h" // For prefetch if needed, though likely not critical here
#include "types.h"

namespace Stockfish {

namespace KpkBitbase {

namespace { // Anonymous namespace for internal details

  // --- Constants and Data Structures (from old implementation) ---

  // There are 24 possible pawn squares: files A to D and ranks from 2 to 7.
  // Positions with the pawn on files E to H will be mirrored before probing.
  constexpr unsigned MAX_INDEX = 2 * 24 * 64 * 64; // stm * psq * wksq * bksq = 196608

  // Using std::bitset for the bitbase storage. Consider alignment if performance critical.
  // alignas(64) // Optional alignment
  std::bitset<MAX_INDEX> KPKBitbase;

  // A KPK bitbase index is an integer in [0, MAX_INDEX) range
  //
  // Information is mapped in a way that minimizes the number of iterations:
  //
  // bit  0- 5: white king square (from SQ_A1 to SQ_H8)
  // bit  6-11: black king square (from SQ_A1 to SQ_H8)
  // bit    12: side to move (WHITE or BLACK)
  // bit 13-14: white pawn file (from FILE_A to FILE_D)
  // bit 15-17: white pawn RANK_7 - rank (from RANK_7 - RANK_7 to RANK_7 - RANK_2)
  //
  // IMPORTANT: This index function assumes the pawn 'psq' belongs to WHITE.
  unsigned index(Color stm, Square bksq, Square wksq, Square psq) {
    assert(file_of(psq) <= FILE_D); // Mirrored pawn square must be A-D file
    assert(rank_of(psq) >= RANK_2 && rank_of(psq) <= RANK_7); // Pawn rank must be valid

    // Use static_cast for clarity converting enums/scoped enums to int/unsigned
    return static_cast<unsigned>(wksq)
         | (static_cast<unsigned>(bksq) << 6)
         | (static_cast<unsigned>(stm) << 12)
         | (static_cast<unsigned>(file_of(psq)) << 13)
         | ((static_cast<unsigned>(RANK_7) - static_cast<unsigned>(rank_of(psq))) << 15);
  }

  enum Result {
    INVALID = 0,
    UNKNOWN = 1,
    DRAW    = 2, // Includes loss for white (stalemate, pawn capture)
    WIN     = 4  // Win for white
  };

  // Overload operator|= for convenient flag setting
  Result& operator|=(Result& r, Result v) { return r = static_cast<Result>(r | v); }

  // Helper struct for initialization (from old implementation)
  struct KPKPosition {
    KPKPosition() = default;
    explicit KPKPosition(unsigned idx); // Constructor to decode index
    operator Result() const { return result; } // Conversion to Result
    Result classify(const std::vector<KPKPosition>& db); // Retrograde analysis step

    Color stm;
    Square ksq[COLOR_NB], psq; // ksq[WHITE], ksq[BLACK], psq (White's pawn)
    Result result;
  };


  // --- Initialization Logic (adapted from old implementation) ---

  // Constructor: Decodes index to position details and sets initial known results
  KPKPosition::KPKPosition(unsigned idx) {

    ksq[WHITE] = static_cast<Square>((idx >>  0) & 0x3F);
    ksq[BLACK] = static_cast<Square>((idx >>  6) & 0x3F);
    stm        = static_cast<Color> ((idx >> 12) & 0x01);
    psq        = make_square(static_cast<File>((idx >> 13) & 0x3),
                             static_cast<Rank>(RANK_7 - ((idx >> 15) & 0x7))); // Rank 7 - offset

    result = UNKNOWN; // Default state

    // Basic validity checks
    if (!is_ok(ksq[WHITE]) || !is_ok(ksq[BLACK]) || !is_ok(psq) // Check squares are valid
        || ksq[WHITE] == ksq[BLACK] // Kings too close is handled by distance check below
        || ksq[WHITE] == psq
        || ksq[BLACK] == psq
        || rank_of(psq) < RANK_2 || rank_of(psq) > RANK_7 // Pawn on invalid rank
       )
    {
        result = INVALID;
        return;
    }

    // Invalid if kings are adjacent or overlapping
    if (distance(ksq[WHITE], ksq[BLACK]) <= 1)
    {
        result = INVALID;
        return;
    }

    // Invalid if the king of the side *not* to move is attacked by the pawn
    // This represents a state that shouldn't be reached if the previous move was legal.
    if (pawn_attacks_bb<WHITE>(psq) & ksq[BLACK]) // White pawn attacks Black king
    {
        if (stm == WHITE) // If it's White's turn, Black King must not be in check by pawn initially
        {
             result = INVALID;
             return;
        }
        // If it's Black's turn and BK is attacked by pawn, it's a WIN for White immediately (already handled?)
        // The original code checks this only for stm == WHITE. Let's stick to that for initial re-implementation.
        // It implies the *previous* state was illegal if BK is attacked and it's now White's turn.
    }

    // --- Initial Known Results ---

    // Win if the pawn can promote safely next move
    Square promotion_sq = psq + pawn_push(WHITE); // Assuming WHITE pawn push direction
    if (stm == WHITE && rank_of(psq) == RANK_7)
    {
        // Promotion square must be valid and different from White King's position
        if (is_ok(promotion_sq) && ksq[WHITE] != promotion_sq)
        {
            // Win if Black King cannot capture the promoting pawn or the square,
            // OR if White King defends the promotion square.
            if (distance(ksq[BLACK], promotion_sq) > 1 || distance(ksq[WHITE], promotion_sq) <= 1)
            {
                result = WIN;
                return;
            }
        }
    }

    // Draw if it's Black to move and it's stalemate or Black king captures the pawn
    if (stm == BLACK)
{
    Bitboard bk_moves = attacks_bb<KING>(ksq[BLACK]);
    // Squares occupied by White pieces (only King and Pawn in KPK)
    Bitboard white_occupied = square_bb(ksq[WHITE]) | square_bb(psq);
    // Squares attacked by White or occupied by White
    Bitboard unsafe_for_black = attacks_bb<KING>(ksq[WHITE]) | pawn_attacks_bb<WHITE>(psq) | white_occupied;
    Bitboard safe_sq = ~unsafe_for_black; // Squares not attacked by White King or Pawn, or occupied by White

    // Stalemate check: No legal moves for Black King
    if (!(bk_moves & safe_sq))
        {
            result = DRAW;
            return;
        }

        // Black King can capture the pawn (and the square is safe from White King)
        if ((attacks_bb<KING>(ksq[BLACK]) & psq) && !(attacks_bb<KING>(ksq[WHITE]) & psq))
        {
             result = DRAW;
             return;
        }
    }
  }


  // Classification function for retrograde analysis
  Result KPKPosition::classify(const std::vector<KPKPosition>& db) {

    // If already classified, return
    if (result != UNKNOWN)
        return result;

    const Result TargetResult = (stm == WHITE ? WIN : DRAW); // What the side to move wants
    const Result FallbackResult = (stm == WHITE ? DRAW : WIN); // What happens if no target move found

    Result currentClassification = INVALID; // Start assuming no legal moves lead anywhere good
    bool hasUnknown = false; // Track if any successor state is still UNKNOWN

    // 1. Generate King Moves
Bitboard b = attacks_bb<KING>(ksq[stm]);
while (b) {
    Square to_sq = pop_lsb(b);

        // Check basic legality (destination square valid and not occupied by friendly piece)
        if (!is_ok(to_sq) || to_sq == ksq[~stm] || to_sq == psq)
             continue; // Invalid destination

        // Check if king moves into check from the other king
        if (distance(to_sq, ksq[~stm]) <= 1)
            continue;

        // Check if king moves into check from the pawn (only if stm is BLACK moving)
        if (stm == BLACK && (pawn_attacks_bb<WHITE>(psq) & to_sq))
            continue;

        // Generate index for the resulting position
        unsigned next_idx;
        if (stm == WHITE)
            next_idx = index(BLACK, ksq[BLACK], to_sq, psq);
        else // stm == BLACK
            next_idx = index(WHITE, to_sq, ksq[WHITE], psq);

        // Check the result from the database
        Result next_result = db[next_idx];

        // If any move leads to the target result, the current position is classified as such
        if (next_result == TargetResult) {
            result = TargetResult;
            return result;
        }

        // Accumulate results
        if (next_result != INVALID) {
             currentClassification |= next_result;
             if (next_result == UNKNOWN) hasUnknown = true;
        }
    }

    // 2. Generate Pawn Moves (only if stm == WHITE)
    if (stm == WHITE) {
        // a. Single Push
        Square push1_sq = psq + pawn_push(WHITE);
        if (rank_of(psq) < RANK_7 && is_ok(push1_sq) && push1_sq != ksq[WHITE] && push1_sq != ksq[BLACK])
        {
            unsigned next_idx = index(BLACK, ksq[BLACK], ksq[WHITE], push1_sq);
            Result next_result = db[next_idx];
            if (next_result == TargetResult) { result = TargetResult; return result; }
            if (next_result != INVALID) {
                 currentClassification |= next_result;
                 if (next_result == UNKNOWN) hasUnknown = true;
            }
        }

        // b. Double Push
        Square push2_sq = push1_sq + pawn_push(WHITE);
        if (rank_of(psq) == RANK_2 && is_ok(push1_sq) && is_ok(push2_sq)
            && push1_sq != ksq[WHITE] && push1_sq != ksq[BLACK] // Path clear
            && push2_sq != ksq[WHITE] && push2_sq != ksq[BLACK]) // Destination clear
        {
            unsigned next_idx = index(BLACK, ksq[BLACK], ksq[WHITE], push2_sq);
            Result next_result = db[next_idx];
            if (next_result == TargetResult) { result = TargetResult; return result; }
            if (next_result != INVALID) {
                 currentClassification |= next_result;
                 if (next_result == UNKNOWN) hasUnknown = true;
            }
        }
    }

    // 3. Final Classification
    if (hasUnknown || currentClassification == INVALID) // If any successor is UNKNOWN or no valid moves found
    {
        result = UNKNOWN; // Cannot classify yet or stalemate/invalid state reached by all moves
                          // Note: INVALID should ideally not happen if initial INVALID states are correct
                          // Stalemate for WHITE leads to DRAW. Stalemate for BLACK handled in constructor.
        if (currentClassification == INVALID && stm == WHITE) result = DRAW; // White has no legal moves -> Draw
    }
    else if (currentClassification & TargetResult) // Should have been caught earlier, but check again
    {
        result = TargetResult;
    }
     else if (currentClassification & UNKNOWN) // Should have been caught by hasUnknown
    {
         result = UNKNOWN;
    }
    else // All reachable states lead to the opponent's desired outcome or are draws (for white)
    {
        result = FallbackResult;
    }

    return result;
  }

} // anonymous namespace


// --- Public Interface Implementation ---

void init() {

  std::vector<KPKPosition> db(MAX_INDEX);
  unsigned idx = 0; // Keep idx declared outside loop if needed later
  bool changed = true; // Use bool instead of unsigned for change tracking

  // Initialize db with known win / draw / invalid positions
  for (idx = 0; idx < MAX_INDEX; ++idx) {
      db[idx] = KPKPosition(idx);
  }

  // Iterate using retrograde analysis until no position changes state
  int iterations = 0; // Optional: track iterations
  while (changed) {
      changed = false;
      for (idx = 0; idx < MAX_INDEX; ++idx) {
          // Only try to classify UNKNOWN positions
          if (db[idx] == UNKNOWN) {
              Result oldResult = db[idx].result; // Store previous result (should be UNKNOWN)
              Result newResult = db[idx].classify(db); // Attempt to classify
              if (newResult != oldResult && newResult != UNKNOWN) { // Check if it got classified
                  changed = true; // A classification was made, need another pass
              }
          }
      }
      iterations++;
      // std::cout << "Iteration " << iterations << " completed." << std::endl; // Debug output
  }
  // std::cout << "KPK Bitbase generation took " << iterations << " iterations." << std::endl;

  // Fill the final KPKBitbase with the decisive WIN results
  KPKBitbase.reset(); // Clear the bitset first
  for (idx = 0; idx < MAX_INDEX; ++idx) {
      if (db[idx] == WIN) { // Only set bits for White wins
          KPKBitbase.set(idx);
      }
  }
}

bool probe(Square wksq, Square wpsq, Square bksq, Color stm) {

    assert(is_ok(wksq) && is_ok(wpsq) && is_ok(bksq));
    assert(rank_of(wpsq) >= RANK_2 && rank_of(wpsq) <= RANK_7); // Pawn must be on valid rank

    Square   probe_wpsq = wpsq;
    Square   probe_wksq = wksq;
    Square   probe_bksq = bksq;
    File     pawn_file = file_of(wpsq);

    // Mirror if pawn is on files E-H
    if (pawn_file > FILE_D) {
        probe_wpsq = flip_file(wpsq);
        probe_wksq = flip_file(wksq);
        probe_bksq = flip_file(bksq);
    }

    // Calculate the index using the (potentially mirrored) squares
    unsigned idx = index(stm, probe_bksq, probe_wksq, probe_wpsq);

    // Check the bitset
    return KPKBitbase.test(idx);
}

} // namespace KpkBitbase

} // namespace Stockfish
