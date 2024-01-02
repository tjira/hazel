#pragma once

#include "system.h"

inline std::string ORCA = R"(
! METHOD BASIS

*xyz CHARGE MULTI
)";

class Orca {
    struct Results {
        double E; Matrix G;
    };
public:
    // constructor nad input creator
    Orca(const System& system, const std::string& method);

    // setters
    void gradient(double step); void hessian(double step);

    // runner and input getter
    std::string getInput() const {return input;} Results run() const;

private:
    std::string directory, input, method; System system;
};
