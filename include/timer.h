#pragma once

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

struct Timer {
    // Type definitions
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> timepoint;
    typedef std::chrono::milliseconds millis;

    // Getters
    static long elapsed(timepoint start);
    static timepoint now();

    // Utilities
    static std::string format(long ms);
};
