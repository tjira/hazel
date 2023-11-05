#pragma once

#include <fstream>
#include <iomanip>
#include <vector>

class Plot {
public:
    Plot(const std::string& name) : name(name) {}
    void plot(const std::vector<std::vector<double>>& data) const;

private:
    std::string name;
};
