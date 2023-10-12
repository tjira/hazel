#pragma once

#include <iomanip>
#include <sstream>
#include <chrono>
#include <string>

namespace Timer {
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> Timepoint;
    typedef std::chrono::milliseconds Millis;

    // getters
    long Elapsed(Timepoint start);
    std::string Local();
    Timepoint Now();

    // utilities
    std::string Format(long ms);
};
