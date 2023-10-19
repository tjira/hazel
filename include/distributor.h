#pragma once

#include <filesystem>

#include "optimizer.h"
#include "utility.h"
#include "lambda.h"
#include "printer.h"
#include "parser.h"
#include "qd.h"
#include "md.h"

#ifndef CXXFLAGS
#define CXXFLAGS "---"
#endif

#ifndef DATADIR
#define DATADIR ""
#endif

class Distributor {
public:
    Distributor(int argc, char** argv);
    ~Distributor(); void run();

private:
    // calculate integrals
    void integrals();

    // RHF distribution
    void rhfrun(); void rhff(const HF::ResultsRestricted& rhfres);
    void rhfg(const HF::ResultsRestricted& rhfres); void rhfo();

    // RHF distribution
    void uhfrun(); void uhff(const HF::ResultsUnrestricted& uhfres);
    void uhfg(const HF::ResultsUnrestricted& uhfres); void uhfo();

    // RCI distribution
    void rcirun(const HF::ResultsRestricted& rhfres); void rcif(const HF::ResultsRestricted& rhfres);
    void rcig(const HF::ResultsRestricted& rhfres); void rcio();

    // RMP2 distribution
    void rmp2run(const HF::ResultsRestricted& rhfres); void rmp2f(const HF::ResultsRestricted& rhfres);
    void rmp2g(const HF::ResultsRestricted& rhfres); void rmp2o();

    // SCAN, MD and QD distribution
    void scan(); void dynamics(); void qdyn();

private:
    // argument parsers, timer and system
    Parser parser; System system;
    Timer::Timepoint start;
};
