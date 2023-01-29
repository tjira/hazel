#include "../include/printer.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/numeric.hpp>

#define STRINGIFY(X) STRING(X)
#define STRING(X) #X

void Printer::printElapsed(int ms) {
    std::cout << "╔══════════════════════════╗\n";
    std::cout << boost::format("║ TOTAL TIME: %s ║") % Timer::format(ms) << "\n";
    std::cout << "╚══════════════════════════╝" << std::endl;
}

void Printer::printIteration(const HartreeFockResult& res, const HartreeFockOptions& opt) {
    boost::format line("║ %4i │ %22.14f │ %.2e │ %.2e │ %6i ║");
    if (res.i == 1) {
        std::cout << "╔════════════════════════════════════════════════════════════════════╗\n";
        std::cout << "║                         STARTING SCF CYCLE                         ║\n";
        std::cout << "╠══════╤════════════════════════╤══════════╤══════════╤══════════════╣\n";
        std::cout << "║  ##  │         ENERGY         │    dE    │    dD    │     TIME     ║\n";
        std::cout << "╟──────┼────────────────────────┼──────────┼──────────┼──────────────╢\n";
    }
    std::cout << line % res.i % res.E % res.dE % res.dD % Timer::format(res.times.iters.at(res.i - 1)) << std::endl;
    if ((res.dE < opt.thresh && res.dD < opt.thresh) || res.i == opt.maxiter) {
        std::cout << "╚══════╧════════════════════════╧══════════╧══════════╧══════════════╝" << std::endl;
    }
}

void Printer::printMethod(const HartreeFockOptions& opt) {
    boost::format formatter("║ MAXITER: %4i, THRESH: %.2e, DIIS: [ENABLED: %1i, START: %2i, KEEP: %2i, DAMP: %4.2f] ║");
    std::cout << "╔══════════════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                     HARTREE-FOCK                                     ║\n";
    std::cout << "╠══════════════════════════════════════════════════════════════════════════════════════╣\n";
    std::cout << formatter % opt.maxiter % opt.thresh % opt.diis.enabled % opt.diis.start % opt.diis.keep % opt.diis.damp << "\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════════════╝" << std::endl;
}

void Printer::printResult(const HartreeFockResult& res) {
    boost::format orbital("║ %4i │ %3.1f │ %22.14f │ %22.14f ║");
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                       ORBITAL ENERGIES                       ║\n";
    std::cout << "╠══════╤═════╤════════════════════════╤════════════════════════╣\n";
    std::cout << "║  ##  │ OCC │         E (Eh)         │         E (eV)         ║\n";
    std::cout << "╟──────┼─────┼────────────────────────┼────────────────────────╢\n";
    for (int i = 0; i < res.Eo.rows(); i++) {
        std::cout << orbital % i % ((i + 1) <= res.nocc ? 2.0 : 0.0) % res.Eo(i) % (res.Eo(i) * EH2EV) << "\n";
    }
    std::cout << "╚══════╧═════╧════════════════════════╧════════════════════════╝\n";
    std::cout << "╔══════════════════════════════════════════════════════╗\n";
    std::cout << boost::format("║ FINAL SINGLE POINT ENERGY: %22.14f Eh ║") % res.E << "\n";
    std::cout << "╚══════════════════════════════════════════════════════╝" << std::endl;
}

void Printer::printInitialTimings(const HartreeFockResult& res) {
    std::cout << "╔═══════════════════════════════════════════════════════════╗\n";
    std::cout << "║          ONE-ELECTRON INTEGRAL COMPUTATION TIMES          ║\n";
    std::cout << "╠══════════════╤══════════════╤══════════════╤══════════════╣\n";
    std::cout << "║   T MATRIX   │   V MATRIX   │   S MATRIX   │   COMBINED   ║\n";
    std::cout << "╟──────────────┼──────────────┼──────────────┼──────────────╢\n";
    std::cout << boost::format("║ %s │ %s │ %s │ %s ║") % Timer::format(res.times.ints.at("T")) % Timer::format(res.times.ints.at("V")) %
    Timer::format(res.times.ints.at("S")) % Timer::format(boost::accumulate(res.times.ints | boost::adaptors::map_values, 0)) << "\n";  
    std::cout << "╚══════════════╧══════════════╧══════════════╧══════════════╝\n";
    std::cout << "╔═════════════════════════════╗\n";
    std::cout << "║         GUESS TIMES         ║\n";
    std::cout << "╠══════════════╤══════════════╣\n";
    std::cout << "║   D MATRIX   │   GUESS E0   ║\n";
    std::cout << "╟──────────────┼──────────────╢\n";
    std::cout << boost::format("║ %s │ %s ║") % Timer::format(res.times.guess.at(0)) % Timer::format(res.times.guess.at(1)) << "\n";
    std::cout << "╚══════════════╧══════════════╝" << std::endl;
}

void Printer::printTitle() {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time(); auto tm = to_tm(now);
    std::cout << "╔══════════════════════════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                      QUANTUM HAZEL SOFTWARE                                      ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════════════════════════╝\n";
    std::cout << "╔═════════════════════════════════╗\n";
    std::cout << std::put_time(&tm, "║ TIMESTAMP: %d-%b-%Y %H:%M:%S ║\n");
    std::cout << "╚═════════════════════════════════╝" << std::endl;
    /* std::cout << STRINGIFY(GPPFLAGS) << std::endl; */
}
