#pragma once

#include "forward.h"
#include <boost/format.hpp>
#include <iomanip>
#include <iostream>

class Printer {
public:
    template <class R, class O>
    static void printIteration(const R& res, const O& opt);
    template <class O>
    static void printMethod(const O& opt);
    template <class R>
    static void printResult(const R& res);
    static void printTitle();

private:
    static std::string repeat(const std::string& str, int i);
    static int w1, w2, w3;
};

template <class R, class O>
void Printer::printIteration(const R& res, const O& opt) {
    boost::format line("║ %4i │ %22.14f │ %.2e │ %.2e │ %6i ║");
    if (res.i == 1) {
        std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
        std::cout << "║                      STARTING SCF CYCLE                      ║\n";
        std::cout << "╠══════╤════════════════════════╤══════════╤══════════╤════════╣\n";
        std::cout << "║  ##  │         energy         │    dE    │    dD    │  time  ║\n";
        std::cout << "╟──────┼────────────────────────┼──────────┼──────────┼────────╢\n";
    }
    std::cout << line % res.i % res.E % res.dE % res.dD % res.times.iters.at(res.i - 1) << std::endl;
    if ((res.dE < opt.thresh && res.dD < opt.thresh) || res.i == opt.maxiter) {
        std::cout << "╚══════╧════════════════════════╧══════════╧══════════╧════════╝" << std::endl;
    }
}

template <class O>
void Printer::printMethod(const O& opt) {
    if (std::is_same<O, HartreeFockOptions>()) {
        std::stringstream optss;
        optss<< std::scientific << std::setprecision(0) << "║ MAXITER=" << opt.maxiter << " THRESH=" << opt.thresh;
        optss << std::fixed << std::setprecision(3) << " DAMP=" << opt.damp;
        std::cout << "╔" << repeat("═", w2 - 2) << "╗\n";
        std::cout << "║" << std::string((w2 - 12) / 2 - 1, ' ') << "HARTREE-FOCK" << std::string((w2 - 12) / 2 - 1, ' ') << "║\n";
        std::cout << "╠" << repeat("═", w2 - 2) << "╣\n";
        std::cout << optss.str() << repeat(" ", w2 - optss.str().size() + 1) << "║\n";
        std::cout << "╚" << repeat("═", w2 - 2) << "╝" << std::endl;
    }
}

template <class R>
void Printer::printResult(const R& res) {
    boost::format final("║ FINAL SINGLE POINT ENERGY: %22.14f Eh ║");
    boost::format orbital("║ %4i │ %3.1f │ %22.14f │ %22.14f ║");
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                       ORBITAL ENERGIES                       ║\n";
    std::cout << "╠══════╤═════╤════════════════════════╤════════════════════════╣\n";
    std::cout << "║  ##  │ occ │         E (Eh)         │         E (eV)         ║\n";
    std::cout << "╟──────┼─────┼────────────────────────┼────────────────────────╢\n";
    for (int i = 0; i < res.Eo.rows(); i++) {
        double occ = ((i + 1) <= res.nocc ? 2.0 : 0.0);
        std::cout << orbital % i % occ % res.Eo(i) % (res.Eo(i) * EH2EV) << "\n";
    }
    std::cout << "╚══════╧═════╧════════════════════════╧════════════════════════╝\n";
    std::cout << "╔══════════════════════════════════════════════════════╗\n";
    std::cout << final % res.E << "\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n";
    std::cout << std::flush;
}
