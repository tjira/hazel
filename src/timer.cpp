#include "timer.h"

long Timer::Elapsed(Timepoint start) {
    return std::chrono::duration_cast<Millis>(Now() - start).count();
};

std::string Timer::Format(long ms) {
    long hours = ms / 3600000, mins = ms % 3600000 / 60000;
    long secs = ms % 60000 / 1000; ms = ms % 1000;
    std::stringstream ss; ss << std::setfill('0');
    ss << std::setw(2) << hours <<  ":" << std::setw(2) << mins << ":";
    ss << std::setw(2) << secs << "." << std::setw(3) << ms;
    return ss.str();
}

std::string Timer::Local() {
    auto t = std::time(nullptr); auto tm = *std::localtime(&t);
    std::stringstream ss; ss << std::put_time(&tm, "%a %b %e %T %Y");
    return ss.str();
}

Timer::Timepoint Timer::Now() {
    return std::chrono::high_resolution_clock().now();
};
