#pragma once

#include <Eigen/Eigen>

#define EH2EV 27.211324570273

class HartreeFock;
class Molecule;

struct HartreeFockOptions {
    double damp, thresh;
    int maxiter;
};

struct HartreeFockResult {
    Eigen::MatrixXd T, V, S; double Vnn;
    std::vector<Eigen::MatrixXd> Fs, Ds;
    std::vector<double> Es;
    Eigen::VectorXd Eo;
    int nocc, iters;
    struct {
        std::vector<long> iters;
        long guess, ints;
    } times;
};
