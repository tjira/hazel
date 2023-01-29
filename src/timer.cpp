#include "../include/timer.h"

typedef std::chrono::time_point<std::chrono::high_resolution_clock> timepoint;
typedef std::chrono::milliseconds millis;

long Timer::elapsed(timepoint start) {
    return std::chrono::duration_cast<millis>(now() - start).count();
};

std::string Timer::format(long ms) {
    long hours = ms / 3600000, mins = ms % 3600000 / 60000;
    long secs = ms % 60000 / 1000; ms = ms % 1000;
    std::stringstream ss; ss << std::setfill('0');
    ss << std::setw(2) << hours <<  ":" << std::setw(2) << mins << ":";
    ss << std::setw(2) << secs << "." << std::setw(3) << ms;
    return ss.str();
}

timepoint Timer::now() {
    return std::chrono::high_resolution_clock().now();
};
