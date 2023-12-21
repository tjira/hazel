#include "plot.h"

void Plot::plot(const std::vector<std::vector<double>>& data) const {
    // open the file
    std::ofstream fstream(name + ".gpl");

    // write name and format to the file
    fstream << "set term " << (term == "png" ? "pngcairo" : term) << " "
            << "size " << size.at(0) << ", " << size.at(1) << "; "
            << "set output \"" << name << "." << term << "\""
            << std::endl;

    // write plot method
    fstream << "plot \"-\" with lines lw " << width << " lc \"blue\" notitle" << std::endl;

    // write the data
    for (const std::vector<double>& row : data) {
        fstream << std::fixed << std::setprecision(14) << std::setw(20) << row.at(0) << " " << row.at(1) << std::endl;
    }
}
