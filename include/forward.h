#pragma once

#include <Eigen/Eigen>
#include <chrono>

class HartreeFock;
class Molecule;

struct HartreeFockOptions {
    double damp, thresh;
    int maxiter;
};

struct HartreeFockResult {
    Eigen::MatrixXd T, V, S; double Vnn;
    std::vector<Eigen::MatrixXd> Fs, Ds;
    std::vector<double> Es, DNs;
    struct {
        std::vector<long> iters;
        long guess, ints;
    } times;
};

struct Timer {
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> timepoint;
    typedef std::chrono::milliseconds millis;

    static long elapsed(timepoint start) { return std::chrono::duration_cast<millis>(now() - start).count(); };
    static timepoint now() { return std::chrono::high_resolution_clock().now(); };
};
