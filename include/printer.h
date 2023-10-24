#pragma once

#include "parser.h"
#include "system.h"

class Printer {
public:
    static void Initial(const Parser& program, System& system);
    static void Mat(const std::string& title, const Matrix& A);
    static void Title(const std::string& title);
};
