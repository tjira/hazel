#pragma once

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

namespace Timer {
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> Timepoint;
    typedef std::chrono::milliseconds Millis;

    // getters
    long Elapsed(Timepoint start);
    Timepoint Now();

    // utilities
    std::string Format(long ms);
};
