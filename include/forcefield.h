#pragma once

#include "potential.h"
#include "system.h"
#include <unordered_map>

struct htuple {
    template <class T1, class T2>
    size_t operator()(const std::tuple<T1, T2>& t) const {
        return std::get<0>(t) ^ std::get<1>(t);
    }
};

class ForceField {
public:
    ForceField(PotentialOptions pairOpt, PotentialOptions bondOpt, System system);
    Eigen::Vector3d F(const std::vector<Particle>& particles, int i) const;
    double U(const std::vector<Particle>& particles) const;

private:
    std::unordered_map<std::tuple<int, int>, std::shared_ptr<Potential>, htuple> pair, bond;
};
