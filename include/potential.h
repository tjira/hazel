#pragma once

#include "particle.h"
#include <include/forward.h>

class Potential {
public:
    virtual Eigen::Vector3d F(const Particle& a, const Particle& b) const { return {0, 0, 0}; };
    virtual double U(const Particle& a, const Particle& b) const { return 0; };
};

struct PotentialArray {
    PotentialArray(std::string name, std::vector<PotentialCoefficients> coefs);
    Eigen::Vector3d F(const Particle& a, const Particle& b) const;
    double U(const Particle& a, const Particle& b) const;
    std::vector<std::vector<std::shared_ptr<Potential>>> array;
};

class LennardJones : public Potential {
public:
    LennardJones(double epsilon, double sigma) : Potential(), epsilon(epsilon), sigma(sigma) {}
    Eigen::Vector3d F(const Particle& a, const Particle& b) const;
    double U(const Particle& a, const Particle& b) const;

private:
    double epsilon, sigma;
};
