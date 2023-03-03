#pragma once

#include "potential.h"
#include "system.h"
#include <unordered_map>

struct htuple {
    template <class T1, class T2>
    size_t operator()(const std::tuple<T1, T2>& t) const {
        return std::get<0>(t) ^ std::get<1>(t);
    }
    template <class T1, class T2, class T3>
    size_t operator()(const std::tuple<T1, T2, T3>& t) const {
        return std::get<0>(t) ^ std::get<1>(t) ^ std::get<2>(t);
    }
};

class ForceField {
public:
    ForceField(PotentialOptions pairOpt, PotentialOptions bondOpt, PotentialOptions angleOpt);
    Eigen::Vector3d F(const std::vector<Particle>& particles, int i) const;
    double U(const std::vector<Particle>& particles) const;

private:
    std::unordered_map<std::tuple<int, int, int>, std::shared_ptr<Potential>, htuple> angle;
    std::unordered_map<std::tuple<int, int>, std::shared_ptr<Potential>, htuple> pair, bond;
};
