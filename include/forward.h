#pragma once

#include <Eigen/Eigen>

#define EH2EV 27.211324570273

class HartreeFock;
class Molecule;

struct HartreeFockOptions {
    double damp, thresh;
    int maxiter;
    struct diis {
        int start, keep;
        double damp;
    } diis;
};

struct HartreeFockResult {
    Eigen::MatrixXd T, V, S, F, D;
    double dD, dE, E, Vnn;
    Eigen::VectorXd Eo;
    int nocc, i;
    struct {
        std::vector<long> iters;
        long guess, ints;
    } times;
};
