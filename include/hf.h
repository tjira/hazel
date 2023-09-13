#pragma once

#include "integral.h"
#include "system.h"

class HF {
public:
    struct OptionsRestricted {
        static HF::OptionsRestricted Load(const argparse::ArgumentParser& parser, bool nocoulomb);
        struct {int start, keep;} diis; double thresh; int maxiter; bool nocoulomb;
    };
    struct ResultsRestricted {
        Matrix C, D; Vector eps; double E, Eel, Enuc;
    };
    struct OptionsUnrestricted {
        static HF::OptionsUnrestricted Load(const argparse::ArgumentParser& parser, bool nocoulomb);
        struct {int start, keep;} diis; double thresh; int maxiter; bool nocoulomb;
    };
    struct ResultsUnrestricted {
        Matrix Ca, Cb, Da, Db; double E, Eel, Enuc; Vector epsa, epsb;
    };
public:
    // constructor
    HF(const OptionsUnrestricted& uopt) : uopt(uopt) {}
    HF(const OptionsRestricted& ropt) : ropt(ropt) {}

    // methods
    ResultsUnrestricted uscf(const System& system, Matrix D, bool print = true) const;
    ResultsRestricted rscf(const System& system, Matrix D, bool print = true) const;

private:
    OptionsUnrestricted uopt;
    OptionsRestricted ropt;
};

inline HF::OptionsRestricted HF::OptionsRestricted::Load(const argparse::ArgumentParser& parser, bool nocoulomb) {
    return {{parser.get<std::vector<int>>("-d").at(0), parser.get<std::vector<int>>("-d").at(1)}, parser.get<double>("-t"), parser.get<int>("-i"), nocoulomb};
}

inline HF::OptionsUnrestricted HF::OptionsUnrestricted::Load(const argparse::ArgumentParser& parser, bool nocoulomb) {
    return {{parser.get<std::vector<int>>("-d").at(0), parser.get<std::vector<int>>("-d").at(1)}, parser.get<double>("-t"), parser.get<int>("-i"), nocoulomb};
}
