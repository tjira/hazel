#pragma once

#include "forward.h"
#include <boost/format.hpp>
#include <iomanip>
#include <iostream>

class Printer {
public:
    template <class R>
    static void printIteration(const R& result, bool last);
    template <class O>
    static void printMethod(const O& opt);
    template <class R>
    static void printResult(const R& result);
    static void printTitle();

private:
    static std::string repeat(const std::string& str, int i);
    static int w1, w2, w3;
};

template <class R>
void Printer::printIteration(const R& result, bool last) {
    boost::format line("║ %4i │ %22.14f │ %.2e │ %.2e │ %6i ║");
    if (result.iters == 1) {
        std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
        std::cout << "║                      STARTING SCF CYCLE                      ║\n";
        std::cout << "╠══════╤════════════════════════╤══════════╤══════════╤════════╣\n";
        std::cout << "║  ##  │         energy         │    dE    │    dD    │  time  ║\n";
        std::cout << "╟──────┼────────────────────────┼──────────┼──────────┼────────╢\n";
    }
    double dD = std::abs(result.Ds.at(result.iters).norm() - result.Ds.at(result.iters - 1).norm());
    double dE = std::abs(result.Es.at(result.iters) - result.Es.at(result.iters - 1));
    std::cout << line % result.iters % result.Es.at(result.iters) % dE % dD % result.times.iters.at(result.iters - 1) << std::endl;
    if (last) std::cout << "╚══════╧════════════════════════╧══════════╧══════════╧════════╝" << std::endl;
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
void Printer::printResult(const R& result) {
    boost::format final("║ FINAL SINGLE POINT ENERGY: %22.14f Eh ║");
    boost::format orbital("║ %4i │ %3.1f │ %22.14f │ %22.14f ║");
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                       ORBITAL ENERGIES                       ║\n";
    std::cout << "╠══════╤═════╤════════════════════════╤════════════════════════╣\n";
    std::cout << "║  ##  │ occ │         E (Eh)         │         E (eV)         ║\n";
    std::cout << "╟──────┼─────┼────────────────────────┼────────────────────────╢\n";
    for (int i = 0; i < result.Eo.rows(); i++) {
        double occ = ((i + 1) <= result.nocc ? 2.0 : 0.0);
        std::cout << orbital % i % occ % result.Eo(i) % (result.Eo(i) * EH2EV) << "\n";
    }
    std::cout << "╚══════╧═════╧════════════════════════╧════════════════════════╝\n";
    std::cout << "╔══════════════════════════════════════════════════════╗\n";
    std::cout << final % result.Es.at(result.Es.size() - 1) << "\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n";
    std::cout << std::flush;
}
