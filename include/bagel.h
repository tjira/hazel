#pragma once

#include "system.h"

inline std::string BAGEL = R"(
{ "bagel" : [

{
  "title" : "molecule",
  "df_basis" : "cc-pvdz-jkfit",
  "angstrom" : "true",
  "geometry" : []
}

]}
)";

class Bagel {
public:
    struct Options {
        std::string method;
    };
    struct Results {
        double E; Vector excs, freq; Matrix G;
    };
public:
    // constructors
    Bagel(const System& system, const Options& opt);

    // setters
    void enableGradient(double step); void enableHessian(double step);

    // runner and input getter
    nlohmann::json getInput() const {return input;} Results run() const;
    std::string getFolder() const {return directory;}

    // extractors
    Vector extractFrequencies(const std::string& output) const;
    Matrix extractGradient(const std::string& output) const;
    Vector extractEnergies(const std::string& output) const;
    double extractEnergy(const std::string& output) const;

private:
    Options opt; System system;
    std::string directory;
    nlohmann::json input;
};
