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

MolecularDynamics::MolecularDynamics(const ForceField& field, Options opt) : field(field), opt(opt) {
    // print method specification
    std::cout << "MOLECULAR DYNAMICS" << std::endl;
    std::cout << boost::format("TIMESTEP: %.4f, STEPS: %i") % opt.timestep % opt.steps;
    std::cout << std::endl << std::endl;
};

Result MolecularDynamics::run(std::vector<Particle> particles, std::string output) const {
    Result result; write(std::ofstream(output), particles, 0);
    std::cout << "ITERATIONS\n  ITER   time [fs]        E [Eh]\n";

    double E = field.U(particles);
    for (size_t j = 0; j < particles.size(); j++) {
        E += particles.at(j).ekin();
    }

    std::cout << boost::format("%8i %9.4f %20.14f") % 0 % 0 % E << std::endl; 

    for (int i = 0; i < opt.steps; i++) {
        auto particlesNew = particles;

        for (size_t j = 0; j < particles.size(); j++) {
            Eigen::Vector3d F = field.F(particles, j);
            particlesNew.at(j).a = -F / particles.at(j).mass;
            particlesNew.at(j).v = particles.at(j).v + 0.5 * (particles.at(j).a + particlesNew.at(j).a) * opt.timestep;
            particlesNew.at(j).q = particles.at(j).q + opt.timestep * (particlesNew.at(j).v + 0.5 * particlesNew.at(j).a * opt.timestep);
        }
        particles = particlesNew;

        double E = field.U(particles);
        for (size_t j = 0; j < particles.size(); j++) {
            E += particles.at(j).ekin();
        }

        std::cout << boost::format("%8i %9.4f %20.14f") % (i + 1) % ((i + 1) * opt.timestep) % E << std::endl; 
        write(std::ofstream(output, std::ios_base::app), particles, i);
    }
    std::cout << std::endl;
    return result;
}
