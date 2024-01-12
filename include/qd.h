#pragma once

#include "parser.h"
#include "system.h"

class QD {
public:
    struct Options {
        static QD::Options Load(const Parser& parser);
        std::string potfile; int iters, nstates;
        double dt, thresh; bool imaginary;
    };
    struct Results {
        std::vector<std::vector<CVector>> states;
        Vector energy; CVector r;
    };
public:
    // constructor
    QD(const Options& opt) : opt(opt) {}

    // methods
    Results run(System system, bool print = true) const;

private:
    Options opt;
};

inline QD::Options QD::Options::Load(const Parser& parser) {
    return {parser.get<std::string>("-f"), parser.get<int>("-i"), parser.get<int>("-n"), parser.get<double>("-s"), parser.get<double>("-t"), parser.get<bool>("--imaginary")};
}
