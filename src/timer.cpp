#include "../include/timer.h"

typedef std::chrono::time_point<std::chrono::high_resolution_clock> timepoint;
typedef std::chrono::milliseconds millis;

long Timer::elapsed(timepoint start) {
    return std::chrono::duration_cast<millis>(now() - start).count();
};


timepoint Timer::now() {
    return std::chrono::high_resolution_clock().now();
};
