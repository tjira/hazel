#pragma once

#include "gradient.h"
#include "hessian.h"

class MD {
public:
    struct Options {
        int iters, state; double step; std::string output;
        struct Berendsen {double temp, tau;} berendsen;
        static MD::Options Load(const Parser& parser);
    };
public:
    // constructor
    MD(const Options& opt) : opt(opt) {}

    // methods
    template <typename F> void run(System system, const F& egfunc, bool print = true) const;

private:
    Options opt;
};

inline MD::Options MD::Options::Load(const Parser& parser) {
    return {parser.get<int>("-i"), parser.get<int>("-e"), parser.get<double>("-s"), parser.get<std::string>("-o"), {BOLTZMANN * parser.get<std::vector<double>>("--berendsen").at(0), parser.get<std::vector<double>>("--berendsen").at(1)}};
}
