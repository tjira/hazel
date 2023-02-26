#include "../include/particle.h"

Particle::Particle(Eigen::Vector3d q, std::string symbol) : q(q), v{0, 0, 0}, a{0, 0, 0}, symbol(symbol) {
    mass = ptable.at(symbol).mass;
}
