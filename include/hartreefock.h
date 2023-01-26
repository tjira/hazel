#pragma once

#include "forward.h"
#include "printer.h"

class HartreeFock {
    typedef HartreeFockOptions Options;
    typedef HartreeFockResult Result;
public:
    HartreeFock(Options opt);
    Result scf(const Molecule& molecule) const;

private:
    Eigen::MatrixXd computeDensity(Eigen::MatrixXd S, Eigen::MatrixXd F, int nocc) const;
    bool checkConvergence(const Result& result, int i) const;
    double computeEnergy(Eigen::ArrayXXd H, Eigen::ArrayXXd F, Eigen::ArrayXXd D) const;
    void logIteration(int i, const Result& result) const;
    Options opt;
};


