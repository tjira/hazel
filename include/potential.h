#pragma once

#include "particle.h"

class Potential {
public:
    virtual Eigen::Vector3d F(const Particle& a, const Particle& b) const { return {0, 0, 0}; };
    virtual double U(const Particle& a, const Particle& b) const { return 0; };
};

class LennardJones : public Potential {
public:
    LennardJones(double epsilon, double sigma) : Potential(), epsilon(epsilon), sigma(sigma) {}
    Eigen::Vector3d F(const Particle& a, const Particle& b) const;
    double U(const Particle& a, const Particle& b) const;

private:
    double epsilon, sigma;
};
