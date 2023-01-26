#pragma once

#include "forward.h"
#include <boost/format.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>

class Printer {
public:
    static void printTitle();
    template <class R> static void printIteration(const R& result, int i, bool last);
    template <class O> static void printMethod(const O& opt);
    template <class R> static void printResult(const R& result);

private:
    static std::string repeat(const std::string& str, int i);
    static int w1, w2, w3;
};

template <class R>
void Printer::printIteration(const R& result, int i, bool last) {
    boost::format line("║ %04i │ %020.14f │ %.2e │ %.2e │ %06i ║");
    if (i == 1) {
        std::cout << "╔════════════════════════════════════════════════════════════╗\n";
        std::cout << "║                     STARTING SCF CYCLE                     ║\n";
        std::cout << "╠══════╤══════════════════════╤══════════╤══════════╤════════╣\n";
        std::cout << "║  ##  │        energy        │    dE    │    dD    │  time  ║\n";
        std::cout << "╟──────┼──────────────────────┼──────────┼──────────┼────────╢\n";
    }
    double dE = std::abs(result.Es.at(i - 1) - result.Es.at(i)), dD = std::abs(result.DNs.at(i - 1) - result.DNs.at(i));
    std::cout << line % i % result.Es.at(i) % dE % dD % result.times.iters.at(i - 1) << std::endl;
    if (last) std::cout << "╚════════════════════════════════════════════════════════════╝" << std::endl;
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
    boost::format energy("║ FINAL SINGLE POINT ENERGY: %020.14f Eh ║");
    std::cout << "╔════════════════════════════════════════════════════╗\n";
    std::cout << energy % result.Es.at(result.Es.size() - 1) << "\n";
    std::cout << "╚════════════════════════════════════════════════════╝\n";
    std::cout << std::flush;
}
