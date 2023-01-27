#pragma once

#include <chrono>

struct Timer {
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> timepoint;
    typedef std::chrono::milliseconds millis;
    static long elapsed(timepoint start);
    static timepoint now();
};
