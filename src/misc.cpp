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

#include "misc.h"

#include <array>
#include <atomic>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <sstream>
#include <string_view>
#include <vector>  // Include for std::vector

#include "types.h"

namespace Stockfish {

namespace {

// Version number or dev.
constexpr std::string_view version = "dev";

// Our fancy logging facility. The trick here is to replace cin.rdbuf() and
// cout.rdbuf() with two Tie objects that tie cin and cout to a file stream. We
// can toggle the logging of std::cout and std:cin at runtime whilst preserving
// usual I/O functionality, all without changing a single line of code!
// Idea from http://groups.google.com/group/comp.lang.c++/msg/1d941c0f26ea0d81

struct Tie: public std::streambuf {  // MSVC requires split streambuf for cin and cout

    Tie(std::streambuf* b, std::streambuf* l) :
        buf(b),
        logBuf(l) {}

    int sync() override { return logBuf->pubsync(), buf->pubsync(); }
    int overflow(int c) override { return log(buf->sputc(char(c)), "<< "); }
    int underflow() override { return buf->sgetc(); }
    int uflow() override { return log(buf->sbumpc(), ">> "); }

    std::streambuf *buf, *logBuf;

    int log(int c, const char* prefix) {

        static int last = '\n';  // Single log file

        if (last == '\n')
            logBuf->sputn(prefix, 3);

        return last = logBuf->sputc(char(c));
    }
};

class Logger {

    Logger() :
        in(std::cin.rdbuf(), file.rdbuf()),
        out(std::cout.rdbuf(), file.rdbuf()) {}
    ~Logger() { start(""); }

    std::ofstream file;
    Tie           in, out;

   public:
    static void start(const std::string& fname) {

        static Logger l;

        if (l.file.is_open())
        {
            std::cout.rdbuf(l.out.buf);
            std::cin.rdbuf(l.in.buf);
            l.file.close();
        }

        if (!fname.empty())
        {
            l.file.open(fname, std::ifstream::out);

            if (!l.file.is_open())
            {
                std::cerr << "Unable to open debug log file " << fname << std::endl;
                exit(EXIT_FAILURE);
            }

            std::cin.rdbuf(&l.in);
            std::cout.rdbuf(&l.out);
        }
    }
};

}  // namespace


// Returns the full name of the current Stockfish version.
//
// For local dev compiles we try to append the commit SHA and
// commit date from git. If that fails only the local compilation
// date is set and "nogit" is specified:
//      Stockfish dev-YYYYMMDD-SHA
//      or
//      Stockfish dev-YYYYMMDD-nogit
//
// For releases (non-dev builds) we only include the version number:
//      Stockfish version
std::string engine_version_info() {
    std::stringstream ss;
    ss << "Stockfish " << version << std::setfill('0');

    if constexpr (version == "dev")
    {
        ss << "-";
#ifdef GIT_DATE
        ss << stringify(GIT_DATE);
#else
        constexpr std::string_view months("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec");

        std::string       month, day, year;
        std::stringstream date(__DATE__);  // From compiler, format is "Sep 21 2008"

        date >> month >> day >> year;
        ss << year << std::setw(2) << std::setfill('0') << (1 + months.find(month) / 4)
           << std::setw(2) << std::setfill('0') << day;
#endif

        ss << "-";

#ifdef GIT_SHA
        ss << stringify(GIT_SHA);
#else
        ss << "nogit";
#endif
    }

    return ss.str();
}

std::string engine_info(bool to_uci) {
    return engine_version_info() + (to_uci ? "\nid author " : " by ")
         + "the Stockfish developers (see AUTHORS file)";
}


// Returns a string trying to describe the compiler we use
std::string compiler_info() {

#define make_version_string(major, minor, patch) \
    stringify(major) "." stringify(minor) "." stringify(patch)

    // Predefined macros hell:
    //
    // __GNUC__                Compiler is GCC, Clang or ICX
    // __clang__               Compiler is Clang or ICX
    // __INTEL_LLVM_COMPILER   Compiler is ICX
    // _MSC_VER                Compiler is MSVC
    // _WIN32                  Building on Windows (any)
    // _WIN64                  Building on Windows 64 bit

    std::string compiler = "\nCompiled by                : ";

#if defined(__INTEL_LLVM_COMPILER)
    compiler += "ICX ";
    compiler += stringify(__INTEL_LLVM_COMPILER);
#elif defined(__clang__)
    compiler += "clang++ ";
    compiler += make_version_string(__clang_major__, __clang_minor__, __clang_patchlevel__);
#elif _MSC_VER
    compiler += "MSVC ";
    compiler += "(version ";
    compiler += stringify(_MSC_FULL_VER) "." stringify(_MSC_BUILD);
    compiler += ")";
#elif defined(__e2k__) && defined(__LCC__)
    #define dot_ver2(n) \
        compiler += char('.'); \
        compiler += char('0' + (n) / 10); \
        compiler += char('0' + (n) % 10);

    compiler += "MCST LCC ";
    compiler += "(version ";
    compiler += std::to_string(__LCC__ / 100);
    dot_ver2(__LCC__ % 100) dot_ver2(__LCC_MINOR__) compiler += ")";
#elif __GNUC__
    compiler += "g++ (GNUC) ";
    compiler += make_version_string(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
    compiler += "Unknown compiler ";
    compiler += "(unknown version)";
#endif

#if defined(__APPLE__)
    compiler += " on Apple";
#elif defined(__CYGWIN__)
    compiler += " on Cygwin";
#elif defined(__MINGW64__)
    compiler += " on MinGW64";
#elif defined(__MINGW32__)
    compiler += " on MinGW32";
#elif defined(__ANDROID__)
    compiler += " on Android";
#elif defined(__linux__)
    compiler += " on Linux";
#elif defined(_WIN64)
    compiler += " on Microsoft Windows 64-bit";
#elif defined(_WIN32)
    compiler += " on Microsoft Windows 32-bit";
#else
    compiler += " on unknown system";
#endif

    compiler += "\nCompilation architecture   : ";
#if defined(ARCH)
    compiler += stringify(ARCH);
#else
    compiler += "(undefined architecture)";
#endif

    compiler += "\nCompilation settings       : ";
    compiler += (Is64Bit ? "64bit" : "32bit");
#if defined(USE_VNNI)
    compiler += " VNNI";
#endif
#if defined(USE_AVX512)
    compiler += " AVX512";
#endif
    compiler += (HasPext ? " BMI2" : "");
#if defined(USE_AVX2)
    compiler += " AVX2";
#endif
#if defined(USE_SSE41)
    compiler += " SSE41";
#endif
#if defined(USE_SSSE3)
    compiler += " SSSE3";
#endif
#if defined(USE_SSE2)
    compiler += " SSE2";
#endif
    compiler += (HasPopCnt ? " POPCNT" : "");
#if defined(USE_NEON_DOTPROD)
    compiler += " NEON_DOTPROD";
#elif defined(USE_NEON)
    compiler += " NEON";
#endif

#if !defined(NDEBUG)
    compiler += " DEBUG";
#endif

    compiler += "\nCompiler __VERSION__ macro : ";
#ifdef __VERSION__
    compiler += __VERSION__;
#else
    compiler += "(undefined macro)";
#endif

    compiler += "\n";

    return compiler;
}


// Debug functions used mainly to collect run-time statistics
// Original: Fixed-size array
// constexpr int MaxDebugSlots = 32;

namespace {

// Original: Fixed-size array
// template<size_t N>
// struct DebugInfo {
//     std::array<std::atomic<int64_t>, N> data = {0};
//
//     [[nodiscard]] constexpr std::atomic<int64_t>& operator[](size_t index) {
//         assert(index < N);
//         return data[index];
//     }
//
//     constexpr DebugInfo& operator=(const DebugInfo& other) {
//         for (size_t i = 0; i < N; i++)
//             data[i].store(other.data[i].load());
//         return *this;
//     }
// };

// Updated: Dynamically sized vector
template<size_t N>
struct DebugInfo {
    std::vector<std::atomic<int64_t>> data;

    DebugInfo(size_t size) : data(size * N) {}  // Initialize with size*N elements

     [[nodiscard]] std::atomic<int64_t>& operator[](size_t index) {
        assert(index < data.size()); //check for the number of elements and not N.
        return data[index];
    }

    DebugInfo& operator=(const DebugInfo& other) {
         // Ensure 'other' has enough elements before copying
        if (other.data.size() >= data.size()) {
            for (size_t i = 0; i < data.size(); i++)
                data[i].store(other.data[i].load());
        }
        return *this;
    }
};


// Original: Fixed size based on MaxDebugSlots
// struct DebugExtremes: public DebugInfo<3> {
//     DebugExtremes() {
//         data[1] = std::numeric_limits<int64_t>::min();
//         data[2] = std::numeric_limits<int64_t>::max();
//     }
// };
struct DebugExtremes : public DebugInfo<3> {
    DebugExtremes(size_t size) : DebugInfo<3>(size) {
        for (size_t i = 0; i < size; ++i) { // Initialize for each slot
            data[i * 3 + 1] = std::numeric_limits<int64_t>::min();
            data[i * 3 + 2] = std::numeric_limits<int64_t>::max();
        }
    }
};
// Original: Fixed-size arrays
// std::array<DebugInfo<2>, MaxDebugSlots>  hit;
// std::array<DebugInfo<2>, MaxDebugSlots>  mean;
// std::array<DebugInfo<3>, MaxDebugSlots>  stdev;
// std::array<DebugInfo<6>, MaxDebugSlots>  correl;
// std::array<DebugExtremes, MaxDebugSlots> extremes;

// Updated:  Use dynamic allocation for debug info.
DebugInfo<2>* hit;
DebugInfo<2>* mean;
DebugInfo<3>* stdev;
DebugInfo<6>* correl;
DebugExtremes* extremes;
size_t numDebugSlots = 0; // Keep track of allocated slots

// Initialization function to allocate memory
void init_debug_info(size_t slots) {
    numDebugSlots = slots;
    hit = new DebugInfo<2>(slots);
    mean = new DebugInfo<2>(slots);
    stdev = new DebugInfo<3>(slots);
    correl = new DebugInfo<6>(slots);
    extremes = new DebugExtremes(slots);
}

// Cleanup function to free allocated memory
void cleanup_debug_info() {
     delete hit;
     delete mean;
     delete stdev;
     delete correl;
     delete extremes;
}

}  // namespace

//No change, using .at() with size check inside functions
void dbg_hit_on(bool cond, int slot) {
    ++(*hit)[slot * 2 + 0]; // Access the correct element for the slot
    if (cond)
        ++(*hit)[slot * 2 + 1];
}

void dbg_mean_of(int64_t value, int slot) {
    ++(*mean)[slot * 2 + 0];
    (*mean)[slot * 2 + 1] += value;
}

void dbg_stdev_of(int64_t value, int slot) {
    ++(*stdev)[slot * 3 + 0];
    (*stdev)[slot * 3 + 1] += value;
    (*stdev)[slot * 3 + 2] += value * value;
}

void dbg_extremes_of(int64_t value, int slot) {
      ++(*extremes)[slot * 3 + 0];

    int64_t current_max = (*extremes)[slot * 3 + 1].load();
    while (current_max < value && !(*extremes)[slot * 3 + 1].compare_exchange_weak(current_max, value))
    {}

    int64_t current_min = (*extremes)[slot * 3 + 2].load();
    while (current_min > value && !(*extremes)[slot * 3 + 2].compare_exchange_weak(current_min, value))
    {}
}

void dbg_correl_of(int64_t value1, int64_t value2, int slot) {
     ++(*correl)[slot * 6 + 0];
    (*correl)[slot * 6 + 1] += value1;
    (*correl)[slot * 6 + 2] += value1 * value1;
    (*correl)[slot * 6 + 3] += value2;
    (*correl)[slot * 6 + 4] += value2 * value2;
    (*correl)[slot * 6 + 5] += value1 * value2;
}

void dbg_print() {

    int64_t n;
    auto    E   = [&n](int64_t x) { return double(x) / n; };
    auto    sqr = [](double x) { return x * x; };

    for (int i = 0; i < numDebugSlots; ++i) // Use numDebugSlots
        if ((n = (*hit)[i*2]))
            std::cerr << "Hit #" << i << ": Total " << n << " Hits " << (*hit)[i*2 + 1]
                      << " Hit Rate (%) " << 100.0 * E((*hit)[i*2 + 1]) << std::endl;

    for (int i = 0; i < numDebugSlots; ++i) // Use numDebugSlots
        if ((n = (*mean)[i*2]))
        {
            std::cerr << "Mean #" << i << ": Total " << n << " Mean " << E((*mean)[i*2 + 1]) << std::endl;
        }

    for (int i = 0; i < numDebugSlots; ++i) // Use numDebugSlots
        if ((n = (*stdev)[i*3]))
        {
            double r = sqrt(E((*stdev)[i*3 + 2]) - sqr(E((*stdev)[i*3 + 1])));
            std::cerr << "Stdev #" << i << ": Total " << n << " Stdev " << r << std::endl;
        }

    for (int i = 0; i < numDebugSlots; ++i) // Use numDebugSlots
        if ((n = (*extremes)[i*3]))
        {
            std::cerr << "Extremity #" << i << ": Total " << n << " Min " << (*extremes)[i*3 + 2]
                      << " Max " << (*extremes)[i*3 + 1] << std::endl;
        }

    for (int i = 0; i < numDebugSlots; ++i) // Use numDebugSlots
        if ((n = (*correl)[i*6]))
        {
            double r = (E((*correl)[i*6 + 5]) - E((*correl)[i*6 + 1]) * E((*correl)[i*6 + 3]))
                     / (sqrt(E((*correl)[i*6 + 2]) - sqr(E((*correl)[i*6 + 1])))
                        * sqrt(E((*correl)[i*6 + 4]) - sqr(E((*correl)[i*6 + 3]))));
            std::cerr << "Correl. #" << i << ": Total " << n << " Coefficient " << r << std::endl;
        }
}
void dbg_clear() {
    //Original Fixed size array
    // hit.fill({});
    // mean.fill({});
    // stdev.fill({});
    // correl.fill({});
    // extremes.fill({});

    // Reset all atomic values to 0 using a loop.  More efficient than recreating the objects.
    for (size_t i = 0; i < hit->data.size(); i++)
        (*hit)[i] = 0;
    for (size_t i = 0; i < mean->data.size(); i++)
         (*mean)[i] = 0;
    for (size_t i = 0; i < stdev->data.size(); i++)
        (*stdev)[i] = 0;
    for (size_t i = 0; i < correl->data.size(); i++)
        (*correl)[i] = 0;
     for (size_t i = 0; i < extremes->data.size(); i++)
        (*extremes)[i] = 0;
}

// Used to serialize access to std::cout
// to avoid multiple threads writing at the same time.
std::ostream& operator<<(std::ostream& os, SyncCout sc) {

    static std::mutex m;

    if (sc == IO_LOCK)
        m.lock();

    if (sc == IO_UNLOCK)
        m.unlock();

    return os;
}

void sync_cout_start() { std::cout << IO_LOCK; }
void sync_cout_end() { std::cout << IO_UNLOCK; }

// Trampoline helper to avoid moving Logger to misc.h
void start_logger(const std::string& fname) { Logger::start(fname); }


#ifdef NO_PREFETCH

void prefetch(const void*) {}

#else

void prefetch(const void* addr) {

    #if defined(_MSC_VER)
    _mm_prefetch((char const*) addr, _MM_HINT_T0);
    #else
    __builtin_prefetch(addr);
    #endif
}

#endif

#ifdef _WIN32
    #include <direct.h>
    #define GETCWD _getcwd
#else
    #include <unistd.h>
    #define GETCWD getcwd
#endif

size_t str_to_size_t(const std::string& s) {
    unsigned long long value = std::stoull(s);
    if (value > std::numeric_limits<size_t>::max())
        std::exit(EXIT_FAILURE);
    return static_cast<size_t>(value);
}

std::optional<std::string> read_file_to_string(const std::string& path) {
    std::ifstream f(path, std::ios_base::binary);
    if (!f)
        return std::nullopt;
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

void remove_whitespace(std::string& s) {
    s.erase(std::remove_if(s.begin(), s.end(), [](char c) { return std::isspace(c); }), s.end());
}

bool is_whitespace(std::string_view s) {
    return std::all_of(s.begin(), s.end(), [](char c) { return std::isspace(c); });
}

std::string CommandLine::get_binary_directory(std::string argv0) {
    std::string pathSeparator;

#ifdef _WIN32
    pathSeparator = "\\";
    #ifdef _MSC_VER
    // Under windows argv[0] may not have the extension. Also _get_pgmptr() had
    // issues in some Windows 10 versions, so check returned values carefully.
    char* pgmptr = nullptr;
    if (!_get_pgmptr(&pgmptr) && pgmptr != nullptr && *pgmptr)
        argv0 = pgmptr;
    #endif
#else
    pathSeparator = "/";
#endif

    // Extract the working directory
    auto workingDirectory = CommandLine::get_working_directory();

    // Extract the binary directory path from argv0
    auto   binaryDirectory = argv0;
    size_t pos             = binaryDirectory.find_last_of("\\/");
    if (pos == std::string::npos)
        binaryDirectory = "." + pathSeparator;
    else
        binaryDirectory.resize(pos + 1);

    // Pattern replacement: "./" at the start of path is replaced by the working directory
    if (binaryDirectory.find("." + pathSeparator) == 0)
        binaryDirectory.replace(0, 1, workingDirectory);

    return binaryDirectory;
}

std::string CommandLine::get_working_directory() {
    std::string workingDirectory = "";
    char        buff[40000];
    char*       cwd = GETCWD(buff, 40000);
    if (cwd)
        workingDirectory = cwd;

    return workingDirectory;
}

}  // namespace Stockfish
