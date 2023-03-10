#include "../include/hartreefock.h"
#include <libint2/diis.h>

static void write(std::ofstream file, System system, int i) {
    file << system.getAtoms().size() << "\n" << "step " << i << "\n";
    for (const auto& atom : system.getAtoms()) {
        file << an2sm.at(atom.atomic_number) << " " << atom.x << " " << atom.y << " " << atom.z << "\n";
    }
}

HartreeFock::MDResult HartreeFock::dynamics(System system, Flags flags) const {
    // create position, velocity and acceleration matrices
    Eigen::MatrixXd a = Eigen::MatrixXd::Zero(system.getAtoms().size(), 3);
    Eigen::MatrixXd q = Eigen::MatrixXd::Zero(system.getAtoms().size(), 3);
    Eigen::MatrixXd v = 0.05 * Eigen::MatrixXd::Random(system.getAtoms().size(), 3);

    // fill the position matrix with the system data
    for (size_t i = 0; i < system.getAtoms().size(); i++) {
        q(i, 0) = system.getAtoms().at(i).x; 
        q(i, 1) = system.getAtoms().at(i).y; 
        q(i, 2) = system.getAtoms().at(i).z; 
    }

    // print the header and write the initial geometry to the file
    Logger::Log(flags.silent, "\nMOLECULAR DYNAMICS ITERATIONS LOOP\n%=8s %=9s %=20s", "ITER", "TIME [fs]", "E [Eh]");
    write(std::ofstream(opt.dyn.output), system, 0);

    // start the MD loop
    for (int i = 0; i < opt.dyn.steps; i++) {

        // start the timer
        auto start = Timer::now();

        // create temporary position, velocity and acc matrices
        Eigen::MatrixXd qp = q, vp = v, ap = a;
    
        // calculate the gradient and energy
        HartreeFock::GDResult grad = gradient(system, { .diis = flags.diis, .silent = true });

        // calculate the acceleration
        for (size_t j = 0; j < system.getAtoms().size(); j++) {
            a(j, 0) = -grad.G(j, 0) / ptable.at(an2sm.at(system.getAtoms().at(j).atomic_number)).mass;
            a(j, 1) = -grad.G(j, 1) / ptable.at(an2sm.at(system.getAtoms().at(j).atomic_number)).mass;
            a(j, 2) = -grad.G(j, 2) / ptable.at(an2sm.at(system.getAtoms().at(j).atomic_number)).mass;
        }

        // calculate the velocity
        v = vp + 0.5 * (ap + a) * opt.dyn.timestep;

        // calculate the difference in position and move the molecule
        Eigen::MatrixXd dq = opt.dyn.timestep * (v + 0.5 * a * opt.dyn.timestep);
        system.move(dq), q = q + dq;
    
        // print the iteration data and write the trajectory
        Logger::Log(flags.silent, "%8i %9.4f %20.14f %s", i + 1, (i + 1) * opt.dyn.timestep, grad.E, Timer::format(Timer::elapsed(start)));
        write(std::ofstream(opt.dyn.output, std::ios_base::app), system, i + 1);
    }

    // return the results
    return {};
};

HartreeFock::GDResult HartreeFock::gradient(System system, Flags flags) const {
    // print the header
    Logger::Log(flags.silent, "\nNUMERICAL GRADIENT ELEMENT EVALUATION LOOP");
    Logger::Log(flags.silent, "%=7s %=17s %=12s", "ELEMENT", "E [Eh/bohr]", "TIME");

    // calculate the energy of the current geometry and initialize gradient matrix
    double E0 = scf(system, { .diis = flags.diis, .silent = true }).E;
    Eigen::MatrixXd G(system.getAtoms().size(), 3);

    // loop over every atom and every coordinate
    #pragma omp parallel for num_threads(opt.engrad.nthread) shared(G)
    for (int i = 0; i < system.getAtoms().size(); i++) {
        for (int j = 0; j < 3; j++) {

            // start the timer
            auto start = Timer::now();

            // create two remporary systems
            System tempSys_minus = system;
            System tempSys_plus = system;

            // offset the temporary systems
            Eigen::MatrixXd dir = Eigen::MatrixXd::Zero(system.getAtoms().size(), 3);
            dir(i, j) = opt.engrad.increment; tempSys_plus.move(dir), tempSys_minus.move(-dir);
            
            // calculate and assign the gradient element
            double E_minus = scf(tempSys_minus, { .diis = flags.diis, .silent = true }).E;
            double E_plus = scf(tempSys_plus, { .diis = flags.diis, .silent = true }).E;
            G(i, j) = (E_plus - E_minus) / opt.engrad.increment / 2;

            // print the element line
            Logger::Log(flags.silent, "(%2i, %1i) %17.14f %s", i + 1, j + 1, G(i, j), Timer::format(Timer::elapsed(start)));
        }
    }

    // print the gradient
    Logger::Log(flags.silent, "\nNUCLEAR GRADIENT");
    Logger::Log(flags.silent, G);

    // return the results
    return { G, E0 };
};

HartreeFock::HFResult HartreeFock::scf(System system, Flags flags) const {
    // calculate the necessary integrals
    Eigen::MatrixXd T = system.integralSingle(libint2::Operator::kinetic);
    Eigen::MatrixXd V = system.integralSingle(libint2::Operator::nuclear);
    Eigen::MatrixXd S = system.integralSingle(libint2::Operator::overlap);
    double Vnn = NUCREP(system); Eigen::MatrixXd H = T + V;

    // print the matrices if requested
    Logger::Log(flags.silent || !opt.print.overlap, "\nOVERLAP MATRIX");
    Logger::Log(flags.silent || !opt.print.overlap, S);
    Logger::Log(flags.silent || !opt.print.kinetic, "\nKINETIC MATRIX");
    Logger::Log(flags.silent || !opt.print.kinetic, T);
    Logger::Log(flags.silent || !opt.print.oneelec, "\nONE-ELECTRON MATRIX");
    Logger::Log(flags.silent || !opt.print.oneelec, V);

    // get the number of occupied orbitals
    int nocc = NELECTRONS(system) / 2;

    // solve the eigenproblem for core hamiltonian
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(H, S);
    Eigen::MatrixXd C = solver.eigenvectors().leftCols(nocc);
    Eigen::VectorXd eps = solver.eigenvalues();

    // compute the guess density matrix and energy
    Eigen::MatrixXd D = 2 * C * C.transpose(); double E = D.cwiseProduct(2 * H).sum() + Vnn;

    // initialize the DIIS algorithm
    libint2::DIIS<Eigen::MatrixXd> diis(opt.diis.start - 1, opt.diis.keep, opt.diis.damp);

    // print the header
    Logger::Log(flags.silent, "\nHARTREE-FOCK SCF CYCLE");
    Logger::Log(flags.silent, "%=4s %=20s %=8s %=8s %=12s", "ITER", "E [Eh]", "dE", "dD", "TIME");

    // start the SCF cycle
    for (int i = 1; i <= opt.maxiter; i++) {

        // start the timer
        auto start = Timer::now();

        // compute the Fock matrix
        Eigen::MatrixXd F = H + 0.5 * system.integralCoulomb(D);

        // compute error and extrapolate the Fock matrix
        Eigen::MatrixXd e = S * D * F - F * D * S;
        if (flags.diis) diis.extrapolate(F, e);

        // save results from previous iterations
        Eigen::MatrixXd Dold = D, Fold = F; double Eold = E;

        // solve the Roothan equations and compute the density matrix from the result
        solver = Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>(F, S);
        C = solver.eigenvectors().leftCols(nocc), eps = solver.eigenvalues();

        // calculate the density matrix and energy
        D = 2 * C * C.transpose(); E = 0.5 * D.cwiseProduct(H + F).sum() + Vnn;

        // calculate the energy change and density norm change
        double dD = std::abs(D.norm() - Dold.norm()), dE = std::abs(E - Eold);

        // print the iteration and restart the timer
        Logger::Log(flags.silent, "%4i %20.14f %.2e %.2e %s",  i, E, dE, dD, Timer::format(Timer::elapsed(start)));

        // check for convergence
        if (dE < opt.thresh && dD < opt.thresh) break;
        else if (i == opt.maxiter) std::cerr << "Algorithm did not converge." << std::endl;
    }
    
    // print the orbital energies if requested
    Logger::Log(flags.silent || !opt.print.orben, "\nORBITAL ENERGIES AND OCCUPATION\n%=4s %=3s %=20s %=22s", "ITER", "OCC", "E [Eh]", "E [eV]");
    for (int i = 0; i < eps.rows(); i++) {
        Logger::Log(flags.silent || !opt.print.orben, "%4i %3.1f %20.14f %22.14f",  i, (i + 1) <= nocc ? 2.0 : 0.0, eps(i), eps(i) * EH2EV);
    }

    // print the orbital coefficients matrix if requested
    Logger::Log(flags.silent || !opt.print.density, "\nORBITAL COEFFICIENTS MATRIX");
    Logger::Log(flags.silent || !opt.print.density, -solver.eigenvectors());

    // print the density matrix if requested
    Logger::Log(flags.silent || !opt.print.density, "\nDENSITY MATRIX");
    Logger::Log(flags.silent || !opt.print.density, D);

    // return the results
    return { C, D, eps, E };
}
