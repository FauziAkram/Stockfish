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

#include "timeman.h"

#include <algorithm>
#include <cmath>

#include "search.h"
#include "uci.h"

namespace Stockfish {
int zz1=50, zz2=100, zz3=120, zz4=330, zz5=44, zz7=200, zz8=680, zz10=1220;
int zz11=84, zz12=335, zz13=30, zz14=480, zz15=360, zz16=300, zz17=270;
int zz18=125 , zz19=1000 , zz20=1120 ;
TUNE(zz1);
TUNE(SetRange(1, 200), zz2);
TUNE(zz3,zz4,zz5,zz7,zz8);
TUNE(SetRange(1, 2400), zz10);
TUNE(zz11,zz12,zz13,zz14,zz15,zz16,zz17,zz18,zz19,zz20);

TimeManagement Time;  // Our global time management object


// Called at the beginning of the search and calculates
// the bounds of time allowed for the current game ply. We currently support:
//      1) x basetime (+ z increment)
//      2) x moves in y seconds (+ z increment)
void TimeManagement::init(Search::LimitsType& limits, Color us, int ply) {

    // If we have no time, no need to initialize TM, except for the start time,
    // which is used by movetime.
    startTime = limits.startTime;
    if (limits.time[us] == 0)
        return;

    TimePoint moveOverhead = TimePoint(Options["Move Overhead"]);
    TimePoint slowMover    = TimePoint(Options["Slow Mover"]);
    TimePoint npmsec       = TimePoint(Options["nodestime"]);

    // optScale is a percentage of available time to use for the current move.
    // maxScale is a multiplier applied to optimumTime.
    double optScale, maxScale;

    // If we have to play in 'nodes as time' mode, then convert from time
    // to nodes, and use resulting values in time management formulas.
    // WARNING: to avoid time losses, the given npmsec (nodes per millisecond)
    // must be much lower than the real engine speed.
    if (npmsec)
    {
        if (!availableNodes)                            // Only once at game start
            availableNodes = npmsec * limits.time[us];  // Time is in msec

        // Convert from milliseconds to nodes
        limits.time[us] = TimePoint(availableNodes);
        limits.inc[us] *= npmsec;
        limits.npmsec = npmsec;
    }

    // Maximum move horizon of 50 moves
    int mtg = limits.movestogo ? std::min(limits.movestogo, zz1) : zz1;

    // Make sure timeLeft is > 0 since we may use it as a divisor
    TimePoint timeLeft = std::max(TimePoint(1), limits.time[us] + limits.inc[us] * (mtg - 1)
                                                  - moveOverhead * (2 + mtg));

    // Use extra time with larger increments
    double optExtra = std::clamp(1.0 + (zz18/10.0) * limits.inc[us] / limits.time[us], (zz19/100.0), (zz20/100.0));

    // Calculate time constants based on current time left.
    double optConstant = std::min((zz12/100000.0) + (zz13/100000.0) * std::log10(limits.time[us] / 1000.0), (zz14/100000.0));
    double maxConstant = std::max((zz15/100.0) + (zz16/100.0) * std::log10(limits.time[us] / 1000.0), (zz17/100.0));

    // A user may scale time usage by setting UCI option "Slow Mover"
    // Default is 100 and changing this value will probably lose elo.
    timeLeft = slowMover * timeLeft / zz2;

    // x basetime (+ z increment)
    // If there is a healthy increment, timeLeft can exceed actual available
    // game time for the current move, so also cap to 20% of available game time.
    if (limits.movestogo == 0)
    {
        optScale = std::min((zz3/10000.0) + std::pow(ply + (zz4/100.0), (zz5/100.0))) * optConstant,
                            (zz7/1000.0) * limits.time[us] / double(timeLeft))
                 * optExtra;
        maxScale = std::min((zz8/100.0), maxConstant + ply / (zz10/100.0));
    }

    // x moves in y seconds (+ z increment)
    else
    {
        optScale = std::min((0.88 + ply / 116.4) / mtg, 0.88 * limits.time[us] / double(timeLeft));
        maxScale = std::min(6.3, 1.5 + 0.11 * mtg);
    }

    // Limit the maximum possible time for this move
    optimumTime = TimePoint(optScale * timeLeft);
    maximumTime =
      TimePoint(std::min((zz11/100.0) * limits.time[us] - moveOverhead, maxScale * optimumTime)) - 10;

    if (Options["Ponder"])
        optimumTime += optimumTime / 4;
}

}  // namespace Stockfish
