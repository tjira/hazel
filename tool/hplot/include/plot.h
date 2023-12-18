#pragma once

#include <fstream>
#include <iomanip>
#include <vector>

class Plot {
public:
    Plot(const std::string& name) : name(name) {}
    void plot(const std::vector<std::vector<double>>& data) const;

private:
    std::vector<int> size = {1024, 768};
    std::string name, term = "png";
    int width = 2;
};
