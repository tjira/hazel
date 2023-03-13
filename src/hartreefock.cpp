#include "../include/hartreefock.h"
#include <libint2/diis.h>

static void write(std::ofstream file, System system, int i) {
    file << system.getAtoms().size() << "\n" << "step " << i << "\n";
    for (const auto& atom : system.getAtoms()) {
        file << an2sm.at(atom.atomic_number) << " " << atom.x << " " << atom.y << " " << atom.z << "\n";
    }
}

HartreeFock::MDResult HartreeFock::dynamics(System system, bool silent) const {
    // print the method header
    Logger::Log(silent, "\nHARTREE-FOCK MOLECULAR DYNAMICS");
    Logger::Log(silent, "TIMESTEP: %.2e, STEPS: %i", opt.dyn.timestep, opt.dyn.steps);

    // create position, velocity and mass matrices
    Mat a = Mat::Zero(system.getSize(), 3), q(system.getSize(), 3), m(system.getSize(), 3);
    Mat v = 0.05 * Mat::Random(system.getAtoms().size(), 3);

    // fill the position and mass matrix
    for(size_t i = 0; i < system.getSize(); i++) {
        m.row(i) = [](double m) { return Vec3(m, m, m); }(ptable.at(an2sm.at(system.getAtom(i).atomic_number)).mass);
        q.row(i) = [](auto a) { return Vec3(a.x, a.y, a.z); }(system.getAtom(i));
    }

    // print the header and write the initial geometry to the file
    Logger::Log(silent, "\nMOLECULAR DYNAMICS ITERATIONS LOOP");
    Logger::Log(silent, "%=8s %=9s %=20s %=12s", "ITER", "TIME [fs]", "E [Eh]", "TIME");
    write(std::ofstream(opt.dyn.output), system, 0);

    // start the MD loop
    for (int i = 0; i < opt.dyn.steps; i++) {

        // start the timer
        auto start = Timer::now();

        // create temporary position, velocity and acc matrices
        Mat qp = q, vp = v, ap = a;
    
        // calculate the gradient and energy
        HartreeFock::GDResult grad = gradient(system, true);

        // perform a Verlet step
        a = -grad.G.array() / m.array();
        v = vp + 0.5 * (ap + a) * opt.dyn.timestep;
        q = qp + opt.dyn.timestep * (v + 0.5 * a * opt.dyn.timestep);

        // update the system
        system = System(system, q);
    
        // print the iteration data and write the trajectory
        Logger::Log(silent, "%8i %9.4f %20.14f %s", i + 1, (i + 1) * opt.dyn.timestep, grad.E, Timer::format(Timer::elapsed(start)));
        write(std::ofstream(opt.dyn.output, std::ios_base::app), system, i + 1);
    }

    // return the results
    return {};
};

HartreeFock::GDResult HartreeFock::gradient(System system, bool silent) const {
    // print the method header
    Logger::Log(silent, "\nHARTREE-FOCK NUMERICAL GRADIENT COMPUTATION");
    Logger::Log(silent, "INCREMENT: %.2e, NTHREAD: %i", opt.engrad.increment, opt.engrad.nthread);

    // print the loop header
    Logger::Log(silent || !opt.engrad.print.iter, "\nHARTREE FOCK NUMERICAL GRADIENT ELEMENT EVALUATION LOOP");
    Logger::Log(silent || !opt.engrad.print.iter, "%=7s %=17s %=12s", "ELEMENT", "E [Eh/bohr]", "TIME");

    // calculate the energy of the current geometry and initialize gradient matrix
    Mat G(system.getAtoms().size(), 3);
    double E0 = scf(system, true).E;

    // loop over every atom and every coordinate
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(opt.engrad.nthread) shared(G)
    #endif
    for (size_t i = 0; i < system.getAtoms().size(); i++) {
        for (size_t j = 0; j < 3; j++) {

            // start the timer
            auto start = Timer::now();

            // create two remporary systems
            System tempSys_minus = system;
            System tempSys_plus = system;

            // offset the temporary systems
            Mat dir = Mat::Zero(system.getAtoms().size(), 3);
            dir(i, j) = opt.engrad.increment; tempSys_plus.move(dir), tempSys_minus.move(-dir);
            
            // calculate and assign the gradient element
            double E_minus = scf(tempSys_minus, true).E;
            double E_plus = scf(tempSys_plus, true).E;
            G(i, j) = (E_plus - E_minus) / opt.engrad.increment / 2;

            // print the element line
            Logger::Log(silent || !opt.engrad.print.iter, "(%2i, %1i) %17.14f %s", i + 1, j + 1, G(i, j), Timer::format(Timer::elapsed(start)));
        }
    }

    // print the gradient
    Logger::Log(silent, "\nNUCLEAR ENERGY GRADIENT");
    Logger::Log(silent, G);

    // return the results
    return { G, E0 };
};

HartreeFock::HFResult HartreeFock::scf(System system, bool silent) const {
    // calculate the necessary integrals
    Mat T = system.integralSingle(libint2::Operator::kinetic);
    Mat V = system.integralSingle(libint2::Operator::nuclear);
    Mat S = system.integralSingle(libint2::Operator::overlap);
    double Vnn = system.getRepulsion(); Mat H = T + V;

    // print the matrices if requested
    Logger::Log(silent || !opt.print.overlap, "\nOVERLAP MATRIX");
    Logger::Log(silent || !opt.print.overlap, S);
    Logger::Log(silent || !opt.print.kinetic, "\nKINETIC MATRIX");
    Logger::Log(silent || !opt.print.kinetic, T);
    Logger::Log(silent || !opt.print.oneelec, "\nONE-ELECTRON MATRIX");
    Logger::Log(silent || !opt.print.oneelec, V);

    // get the number of occupied orbitals
    int nocc = system.getElectrons() / 2;

    // solve the eigenproblem for core hamiltonian
    Eigen::GeneralizedSelfAdjointEigenSolver<Mat> solver(H, S);
    Mat C = solver.eigenvectors().leftCols(nocc);
    Vec eps = solver.eigenvalues();

    // compute the guess density matrix and energy
    Mat D = 2 * C * C.transpose(); double E = D.cwiseProduct(2 * H).sum() + Vnn;

    // initialize the DIIS algorithm
    libint2::DIIS<Mat> diis(opt.diis.start - 1, opt.diis.keep, opt.diis.damp);

    // print the Hartree-Fock header
    Logger::Log(silent, "\nHARTREE-FOCK METHOD");
    if (opt.diis.on) Logger::Log(silent, "MAXITER: %i, THRESH: %.2e, DIIS: [START: %i, KEEP: %i, DAMP: %i]", opt.maxiter, opt.thresh, opt.diis.start, opt.diis.keep, opt.diis.damp);
    else Logger::Log(silent, "MAXITER: %i, THRESH: %.2e, DIIS: OFF", opt.maxiter, opt.thresh);

    // print the SCF header
    Logger::Log(silent || !opt.print.iter, "\nHARTREE-FOCK SCF CYCLE");
    Logger::Log(silent || !opt.print.iter, "%=4s %=20s %=8s %=8s %=12s", "ITER", "E [Eh]", "dE", "dD", "TIME");

    // start the SCF cycle
    for (int i = 1; i <= opt.maxiter; i++) {

        // start the timer
        auto start = Timer::now();

        // compute the Fock matrix and error
        Mat F = H + 0.5 * system.integralCoulomb(D);
        Mat e = S * D * F - F * D * S;

        // compute error and extrapolate the Fock matrix
        if (opt.diis.on) diis.extrapolate(F, e);

        // save results from previous iterations
        Mat Dold = D, Fold = F; double Eold = E;

        // solve the Roothan equations and compute the density matrix from the result
        solver = Eigen::GeneralizedSelfAdjointEigenSolver<Mat>(F, S);
        C = solver.eigenvectors().leftCols(nocc), eps = solver.eigenvalues();

        // calculate the density matrix and energy
        D = 2 * C * C.transpose(); E = 0.5 * D.cwiseProduct(H + F).sum() + Vnn;

        // calculate the energy change and density norm change
        double dD = std::abs(D.norm() - Dold.norm()), dE = std::abs(E - Eold);

        // print the iteration and restart the timer
        Logger::Log(silent || !opt.print.iter, "%4i %20.14f %.2e %.2e %s",  i, E, dE, dD, Timer::format(Timer::elapsed(start)));

        // check for convergence
        if (dE < opt.thresh && dD < opt.thresh) break;
        else if (i == opt.maxiter) std::cerr << "Algorithm did not converge." << std::endl;
    }

    // print the orbital coefficients matrix if requested
    Logger::Log(silent || !opt.print.mos, "\nORBITAL COEFFICIENTS MATRIX");
    Logger::Log(silent || !opt.print.mos, -solver.eigenvectors());

    // print the density matrix if requested
    Logger::Log(silent || !opt.print.density, "\nDENSITY MATRIX");
    Logger::Log(silent || !opt.print.density, D);
    
    // print the orbital energies if requested
    Logger::Log(silent || !opt.print.orben, "\nORBITAL ENERGIES AND OCCUPATION\n%=4s %=3s %=20s %=22s", "ITER", "OCC", "E [Eh]", "E [eV]");
    for (int i = 0; i < eps.rows(); i++) {
        Logger::Log(silent || !opt.print.orben, "%4i %3.1f %20.14f %22.14f",  i, (i + 1) <= nocc ? 2.0 : 0.0, eps(i), eps(i) * EH2EV);
    }

    // print the final energy
    Logger::Log(silent, "\nFINAL SINGLE POINT ENERGY: %.14f Eh", E);

    // return the results
    return { C, D, eps, E };
}
