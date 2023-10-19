#pragma once

#include "parser.h"
#include "system.h"

namespace argparse {
    class ArgumentParser;
};

class Printer {
public:
    static void Initial(const argparse::ArgumentParser& program, System& system);
    static void Mat(const std::string& title, const Matrix& A);
    static void Title(const std::string& title);
};
