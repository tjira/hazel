#pragma once

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

struct Timer {
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> timepoint;
    typedef std::chrono::milliseconds millis;
    static long elapsed(timepoint start);
    static std::string format(long ms);
    static timepoint now();
};
