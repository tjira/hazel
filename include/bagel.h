#pragma once

#include "system.h"

#include <nlohmann/json.hpp>

inline std::string BAGEL = R"(
{ "bagel" : [

{
  "title" : "molecule",
  "df_basis" : "cc-pvdz-jkfit",
  "angstrom" : true,
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
        double E; Vector excs, freq; std::vector<Matrix> Gs;
    };
public:
    // constructors
    Bagel(const System& system, const Options& opt);

    // setters
    void enableGradient(const std::vector<int>& targets = {0}); void enableHessian(double step);

    // runner and input getter
    nlohmann::json getInput() const {return input;} Results run() const;
    std::string getFolder() const {return directory;}

    // extractors
    std::vector<Matrix> extractGradient(const std::string& output) const;
    Vector extractFrequencies(const std::string& output) const;
    Vector extractEnergies(const std::string& output) const;
    double extractEnergy(const std::string& output) const;

private:
    Options opt; System system;
    std::vector<int> targets;
    std::string directory;
    nlohmann::json input;
};
