#pragma once

#include "particle.h"

class Potential {
public:
    virtual Eigen::Vector3d F(const Particle& a, const Particle& b, const Particle& c, int i) const { return {0, 0, 0}; };
    virtual Eigen::Vector3d F(const Particle& a, const Particle& b, int i) const { return {0, 0, 0}; };
    virtual double U(const Particle& a, const Particle& b, const Particle& c) const { return 0; };
    virtual double U(const Particle& a, const Particle& b) const { return 0; };
};

class AngleHarmonic : public Potential {
public:
    AngleHarmonic(double phi0, double k) : Potential(), phi0(phi0), k(k) {}
    Eigen::Vector3d F(const Particle& a, const Particle& b, const Particle& c, int i) const;
    double U(const Particle& a, const Particle& b, const Particle& c) const;

private:
    double phi0, k;
};

class BondHarmonic : public Potential {
public:
    BondHarmonic(double r0, double k) : Potential(), r0(r0), k(k) {}
    Eigen::Vector3d F(const Particle& a, const Particle& b, int i) const;
    double U(const Particle& a, const Particle& b) const;

private:
    double r0, k;
};

class LennardJones : public Potential {
public:
    LennardJones(double epsilon, double sigma) : Potential(), epsilon(epsilon), sigma(sigma) {}
    Eigen::Vector3d F(const Particle& a, const Particle& b, int i) const;
    double U(const Particle& a, const Particle& b) const;

    double epsilon, sigma;
private:
};
