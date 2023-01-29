#pragma once

#include "forward.h"
#include "timer.h"
#include <numeric>
#include <iostream>

namespace Printer {
    void printElapsed(int ms);
    void printIteration(const HartreeFockResult& res, const HartreeFockOptions& opt);
    void printMethod(const HartreeFockOptions& opt);
    void printResult(const HartreeFockResult& res);
    void printInitialTimings(const HartreeFockResult& res);
    void printTitle();
};
