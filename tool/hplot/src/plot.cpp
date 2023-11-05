#include "plot.h"

void Plot::plot(const std::vector<std::vector<double>>& data) const {
    // open the file
    std::ofstream fstream(name + ".gpl");

    // write name and format to the file
    fstream << "set term eps; set output \"" << name << ".eps\"" << std::endl;

    // write plot method
    fstream << "plot \"-\" with lines" << std::endl;

    // write the data
    for (const std::vector<double>& row : data) {
        fstream << std::fixed << std::setprecision(14) << std::setw(20) << row.at(0) << " " << row.at(1) << std::endl;
    }
}
