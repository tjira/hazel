#pragma once

#include <argparse/argparse.hpp>
#include <nlohmann/json.hpp>
#include <Eigen/Eigen>
#include <unordered_map>

#define EH2EV 27.211324570273

using json = nlohmann::json;

class MolecularDynamics;
class HartreeFock;
class Molecule;
class LennardJones;

namespace libint2 {
    inline int nthreads;
}

struct HartreeFockOptions {
    struct DIIS {
        int start, keep;
        bool enabled;
        double damp;
    };
    double thresh; int maxiter; DIIS diis;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HartreeFockOptions::DIIS, start, keep, enabled, damp);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HartreeFockOptions, thresh, maxiter);

struct HartreeFockResult {
    Eigen::MatrixXd T, V, S, F, D, C;
    double dD, dE, E, Vnn;
    Eigen::VectorXd eps;
    int nocc, i;
    struct {
        std::unordered_map<std::string, long> ints;
        std::array<long, 2> guess;
        std::vector<long> iters;
    } times;
};

struct MolecularDynamicsOptions {
    double timestep;
    int steps;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MolecularDynamicsOptions, timestep, steps);

struct MolecularDynamicsResult {

};

struct MullikenResult {
    Eigen::MatrixXd DS;
    Eigen::VectorXd q;
};

struct PotentialCoefficients {
    std::vector<double> coefs;
    std::vector<int> atoms;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PotentialCoefficients, atoms, coefs);

namespace Defaults {
    static json hfopt = R"({
        "name" : "HF",
        "maxiter" : 100,
        "thresh" : 1e-8,
        "diis" : {
            "enabled" : true,
            "start" : 3,
            "keep" : 5,
            "damp" : 0
        }
    })"_json;
    static json mdopt = R"({
        "name" : "MD",
        "output-trajectory" : "trajectory.xyz"
    })"_json;
}
