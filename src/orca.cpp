#include "orca.h"
#include <cstdlib>

Orca::Orca(const System& system, const std::string& method) : input(ORCA), method(method), system(system) {
    #ifndef _WIN32
    // define the directory name and remove the starting new line
    directory = ".orca." + std::to_string(Timer::Now().time_since_epoch().count()); input.erase(0, 1);

    // replace fixed placeholders
    input = std::regex_replace(input, std::regex("BASIS"), Utility::ToUpper(system.basis));
    input = std::regex_replace(input, std::regex("CHARGE"), std::to_string(system.charge));
    input = std::regex_replace(input, std::regex("MULTI"), std::to_string(system.multi));
    input = std::regex_replace(input, std::regex("METHOD"), Utility::ToUpper(method));

    // write the coordinates
    for (int i = 0; i < system.coords.rows(); i++) {
        input += an2sm.at(system.atoms.at(i).atomic_number);
        for (int j = 0; j < 3; j++) {
            input += " " + Utility::ToDblStr(system.coords(i, j));
        } input += "\n";
    }

    // add the star
    input += "*\n";
    #endif
}

void Orca::gradient(double step) {
    if (step) input = std::regex_replace(input, std::regex("\n"), " ENGRAD NUMGRAD\n", std::regex_constants::format_first_only);
    else input = std::regex_replace(input, std::regex("\n"), " ENGRAD\n", std::regex_constants::format_first_only);
}

void Orca::hessian(double step) {
    if (step) input = std::regex_replace(input, std::regex("\n"), " FREQ NUMFREQ\n", std::regex_constants::format_first_only);
    else input = std::regex_replace(input, std::regex("\n"), " FREQ\n", std::regex_constants::format_first_only);
}

Orca::Results Orca::run() const {
    // create directory and open the file
    std::filesystem::create_directory(directory);
    std::ofstream ifile(directory + "/orca.inp");

    // write the input to the opened file
    ifile << input; ifile.close();

    // define the buffer for output and result string with struct
    std::array<char, 128> buffer; std::string result; Results results;

    // create the result matrices
    results.G = Matrix::Zero(system.coords.rows(), 3);

    #ifndef _WIN32
    // execute the command
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(("cd " + directory + " && orca orca.inp | tee orca.out").c_str(), "r"), pclose);

    // check for success
    if (!pipe) throw std::runtime_error("ORCA EXECUTION FAILED");

    // read the output
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += std::string(buffer.data());
    }

    // copy the output
    std::filesystem::copy_file(directory + "/orca.out", "orca.out", std::filesystem::copy_options::overwrite_existing);

    // remove the directory and create thr regex match object
    std::filesystem::remove_all(directory); std::smatch match; std::string searchstr = result, lastmatch;

    // find the energy
    while (std::regex_search(searchstr, match, std::regex("FINAL SINGLE POINT ENERGY(.*?)\n"))) {
        lastmatch = match[1]; searchstr = match.suffix().str();
        try {
            results.E = std::stod(lastmatch); break;
        } catch(...) {};
    }
    
    // create the string stream and line buffer
    std::stringstream lss; lss << result; std::string line;

    // loop over lines
    while (std::getline(lss, line)) {
        if (line == "CARTESIAN GRADIENT") {
            // skip lines
            for (int i = 0; i < 2; i++) std::getline(lss, line);

            // while iterations
            while (std::getline(lss, line) && !line.empty()) {
                // create the column stringstream
                std::stringstream css; css << line;

                // set precision
                css << std::fixed; css.precision(14);

                // exctract data
                int i; std::string str; double x, y, z; css >> i, css >> str, css >> str, css >> x, css >> y, css >> z;

                // set the data
                results.G(i - 1, 0) = x; results.G(i - 1, 1) = y; results.G(i - 1, 2) = z;
            }
        } else if (line == "CARTESIAN GRADIENT (NUMERICAL)") {
            // skip lines
            for (int i = 0; i < 1; i++) std::getline(lss, line);

            // while iterations
            while (std::getline(lss, line) && !line.empty()) {
                // create the column stringstream
                std::stringstream css; css << line;

                // set precision
                css << std::fixed; css.precision(14);

                // exctract data
                int i; std::string str; double x, y, z; css >> i, css >> str, css >> str, css >> x, css >> y, css >> z;

                // set the data
                results.G(i - 1, 0) = x; results.G(i - 1, 1) = y; results.G(i - 1, 2) = z;
            }
        } else if (line == "The final MP2 gradient") {
            // while iterations
            while (std::getline(lss, line) && !line.empty()) {
                // create the column stringstream
                std::stringstream css; css << line;

                // set precision
                css << std::fixed; css.precision(14);

                // exctract data
                std::string str; double x, y, z; css >> str, css >> x, css >> y, css >> z;

                // get the row index
                int i = std::stoi(str.substr(0, str.size() - 1));

                // set the data
                results.G(i, 0) = x; results.G(i, 1) = y; results.G(i, 2) = z;
            }
        }
    }
    #endif

    // return the results
    return results;
}
