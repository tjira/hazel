#pragma once

#include <Eigen/Eigen>
#include <unordered_map>

#define EH2EV 27.211324570273

class MolecularDynamics;
class HartreeFock;
class Molecule;

namespace libint2 {
    inline int nthreads;
}

struct HartreeFockOptions {
    double thresh;
    int maxiter;
    struct diis {
        int start, keep;
        bool enabled;
        double damp;
    } diis;
};

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

struct MolecularDynamicsResult {

};

struct MullikenResult {
    Eigen::MatrixXd DS;
    Eigen::VectorXd q;
};
