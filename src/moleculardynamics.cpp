#include "../include/moleculardynamics.h"

typedef MolecularDynamicsOptions Options;
typedef MolecularDynamicsResult Result;

void write(std::ofstream file, std::vector<Particle> particles, int i) {
    file << particles.size() << "\n" << "step 0\n";
    for (const Particle& particle : particles) {
        Eigen::Vector3d q = particle.getCoords();
        file << particle.getSymbol() << " " << q[0] << " " << q[1] << " " << q[2] << "\n";
    }
}

Result MolecularDynamics::run(std::vector<Particle> particles, std::string output) const {
    Result result; write(std::ofstream(output), particles, 0);

    for (int i = 0; i < opt.steps; i++) {
        auto particlesNew = particles;

        for (size_t j = 0; j < particles.size(); j++) {
            Eigen::Vector3d F = Eigen::Vector3d::Zero();
            for (size_t k = 0; k < particles.size(); k++) {
                if(j != k) F += pot->F(particles.at(j), particles.at(k));
            }
            particlesNew.at(j).a = -F / particles.at(j).mass;
            particlesNew.at(j).v = particles.at(j).v + 0.5 * (particles.at(j).a + particlesNew.at(j).a) * opt.timestep;
            particlesNew.at(j).q = particles.at(j).q + opt.timestep * (particlesNew.at(j).v + 0.5 * particlesNew.at(j).a * opt.timestep);
        }
        particles = particlesNew;

        write(std::ofstream(output, std::ios_base::app), particles, i);
    }
    return result;
}
