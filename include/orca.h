#pragma once

#include "system.h"

inline std::string ORCA = R"(
! METHOD BASIS HCORE NOFROZENCORE

*xyz CHARGE MULTI
)";

class Orca {
public:
    struct Options {
        std::string method;
    };
    struct Results {
        double E; Vector excs, freq; Matrix G;
    };
public:
    // constructors
    Orca(const System& system, const Options& opt);

    // setters
    void enableGradient(double step); void enableHessian(double step);

    // runner and input getter
    std::string getInput() const {return input;} Results run() const;

    // extractors
    Vector extractFrequencies(const std::string& output) const;
    Matrix extractGradient(const std::string& output) const;
    Vector extractEnergies(const std::string& output) const;
    double extractEnergy(const std::string& output) const;

private:
    std::string directory, input;
    Options opt; System system;
};
