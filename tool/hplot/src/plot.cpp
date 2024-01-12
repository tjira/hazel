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
    for (size_t i = 1; i < data.at(0).size(); i++) {
        fstream << (i > 1 ? "     " : "plot ");
        fstream << "\"-\" with lines lw " << width << " lc \"" << colors.at(i - 1) << "\" title \"" << i << "\"";
        fstream << (i < data.at(0).size() - 1 ? ", \\\n" : "\n");
    }

    // write the data
    for (size_t i = 1; i < data.at(0).size(); i++) {
        for (const std::vector<double>& row : data) {
            fstream << std::fixed << std::setprecision(14) << std::setw(20) << row.at(0) << " " << std::setw(20) << row.at(i) << std::endl;
        }
        fstream << "e" << std::endl;
    }
}
