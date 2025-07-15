#include <iostream>

#include "bitboard.h"
#include "misc.h"
#include "position.h"
#include "types.h"
#include "uci.h"
#include "tune.h"
#include "bitbases.h"
#include "psqt.h"
#include "endgame.h"

using namespace Stockfish;

int main(int argc, char* argv[]) {

    std::cout << engine_info() << std::endl;

    Bitboards::init();
    Position::init();
    PSQT::init();
    Bitbases::init();
    Endgames::init();

    UCIEngine uci(argc, argv);

    Tune::init(uci.engine_options());

    uci.loop();

    return 0;
}
