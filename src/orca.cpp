#include "orca.h"
#include "utility.h"

Orca::Orca(const System& system, const Options& opt) : input(ORCA), opt(opt), system(system) {
    // define the directory name and remove the starting new line
    directory = ".orca." + std::to_string(Timer::Now().time_since_epoch().count()); input.erase(0, 1);

    // replace fixed placeholders for charge, multiplicity and basis
    input = std::regex_replace(input, std::regex("BASIS"), Utility::ToUpper(system.originbasis));
    input = std::regex_replace(input, std::regex("CHARGE"), std::to_string(system.charge));
    input = std::regex_replace(input, std::regex("MULTI"), std::to_string(system.multi));

    // replace placeholders for method
    if (Utility::StringContains(opt.method, "casscf")) {
        input = std::regex_replace(input, std::regex("METHOD"), "HF"); std::vector<std::string> casopt; size_t last = 0, next = 0;
        while ((next = opt.method.find("/", last)) != std::string::npos) {
            casopt.push_back(opt.method.substr(last, next - last)); last = next + 1;
        } casopt.push_back(opt.method.substr(last));
        input = std::regex_replace(input, std::regex("\n"), "\n\n%casscf\nnel " + casopt.at(1) + "\nnorb " + casopt.at(2) + "\nmult 1\nnroots " + casopt.at(3) + "\nend\n", std::regex_constants::format_first_only);
    } else {
        input = std::regex_replace(input, std::regex("METHOD"), Utility::ToUpper(opt.method));
    }

    // write the coordinates
    for (int i = 0; i < system.coords.rows(); i++) {
        input += an2sm.at(system.atoms.at(i).atomic_number);
        for (int j = 0; j < 3; j++) {
            input += " " + Utility::ToDblStr(system.coords(i, j));
        } input += "\n";
    }

    // add the star
    input += "*\n";
}

void Orca::enableGradient(double step) {
    if (step) input = std::regex_replace(input, std::regex("\n"), " ENGRAD NUMGRAD\n", std::regex_constants::format_first_only);
    else input = std::regex_replace(input, std::regex("\n"), " ENGRAD\n", std::regex_constants::format_first_only);
}

void Orca::enableHessian(double step) {
    if (step) input = std::regex_replace(input, std::regex("\n"), " FREQ NUMFREQ\n", std::regex_constants::format_first_only);
    else input = std::regex_replace(input, std::regex("\n"), " FREQ\n", std::regex_constants::format_first_only);
}

Orca::Results Orca::run() const {
    // create directory and open the file
    std::filesystem::create_directory(directory);
    std::ofstream ifile(directory + "/orca.inp");

    // write the input to the opened file and define regex match
    ifile << input; ifile.close(); std::smatch match;

    // define the buffer for output and result string with struct
    std::array<char, 128> buffer; std::string output; Results results;

    // execute the command
    #ifdef _WIN32
    std::unique_ptr<FILE, decltype(&pclose)> pipe(_popen(("cd " + directory).c_str(), "r"), _pclose);
    #else
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(("cd " + directory + " && orca orca.inp > >(tee orca.out) 2> /dev/null").c_str(), "r"), pclose);
    #endif

    // check for success
    if (!pipe) throw std::runtime_error("ORCA EXECUTION FAILED");

    // read the output
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        if (Utility::StringContains(std::string(buffer.data()), ": Error")) {
            throw std::runtime_error(std::string(buffer.data()));
        } else if (Utility::StringContains(std::string(buffer.data()), "INPUT ERROR")) {
            throw std::runtime_error(std::string(buffer.data()));
        }
        output += std::string(buffer.data());
    }

    // copy the output
    std::filesystem::copy_file(directory + "/orca.out", "orca.out", std::filesystem::copy_options::overwrite_existing);

    // remove the calculataion directory
    std::filesystem::remove_all(directory);
    
    // return the results
    return {extractEnergy(output), extractEnergies(output), extractFrequencies(output), extractGradient(output)};
}

Vector Orca::extractEnergies(const std::string& output) const {
    // create the string stream and line buffer
    std::stringstream lss; lss << output; std::string line;

    // create the energy vector and some variables
    std::vector<double> excs; std::string str; double energy;

    // loop over lines
    while (std::getline(lss, line)) {
        if (Utility::StringContains(line, "CAS-SCF STATES FOR BLOCK")) {
            // loop over excitation lines
            while (std::getline(lss, line)) {
                if (Utility::StringContains(line, "E=") && !line.empty()) {
                    // create the column stringstream
                    std::stringstream css; css << line;

                    // exctract the correct column energy
                    for (int i = 0; i < 3; i++) {css >> str;} css >> energy;

                    // append the energy
                    excs.push_back(energy);
                }
            } return Eigen::Map<Vector>(excs.data(), excs.size());
        }
    }

    // return the results
    return Eigen::Map<Vector>(excs.data(), excs.size());
}

double Orca::extractEnergy(const std::string& output) const {
    // create the string stream and line buffer
    std::stringstream lss; lss << output; std::string line;

    // loop over lines
    while (std::getline(lss, line)) {
        if (Utility::StringContains(line, "FINAL SINGLE POINT ENERGY")) {
            // create the column stringstream
            std::stringstream css; css << line;

            // set precision
            css << std::fixed; css.precision(14);

            // exctract the correct column energy
            std::string cell; for (int i = 0; i < 5; i++) css >> cell;

            // return the final energy
            if (cell != "Root") return std::stod(cell);
        }
    }

    // return the results
    return 0;
}

Vector Orca::extractFrequencies(const std::string& output) const {
    // create the string stream and line buffer
    std::string line; std::vector<double> f;
    std::stringstream lss; lss << output;

    // create some variables
    std::string str; double freq;

    // loop over lines
    while (std::getline(lss, line)) {
        if (line == "VIBRATIONAL FREQUENCIES") {
            // skip lines
            for (int i = 0; i < 4; i++) std::getline(lss, line);

            // loop over frequency lines
            while (std::getline(lss, line) && !line.empty()) {
                // create the column stringstream
                std::stringstream css; css << line;

                // exctract data
                css >> str, css >> freq;

                // set the data
                if (freq) f.push_back(freq);
            } return Eigen::Map<Vector>(f.data(), f.size());
        }
    }

    // return the frequencies
    return Eigen::Map<Vector>(f.data(), f.size());
}

Matrix Orca::extractGradient(const std::string& output) const {
    // create the string stream and line buffer
    std::stringstream lss; lss << output; std::string line;

    // create the gradient matrix and some variables
    Matrix G = Matrix::Zero(system.coords.rows(), 3);
    int i; double x, y, z; std::string str;

    // loop over lines
    while (std::getline(lss, line)) {
        if (line == "CARTESIAN GRADIENT") {
            // skip lines
            for (int i = 0; i < 2; i++) std::getline(lss, line);

            // loop over gradient lines
            while (std::getline(lss, line) && !line.empty()) {
                // create the column stringstream
                std::stringstream css; css << line;

                // exctract data
                css >> i, css >> str, css >> str, css >> x, css >> y, css >> z;

                // set the data
                G(i - 1, 0) = x; G(i - 1, 1) = y; G(i - 1, 2) = z;
            } return G;
        } else if (line == "CARTESIAN GRADIENT (NUMERICAL)") {
            // skip lines
            for (int i = 0; i < 1; i++) std::getline(lss, line);

            // loop over gradient lines
            while (std::getline(lss, line) && !line.empty()) {
                // create the column stringstream
                std::stringstream css; css << line;

                // exctract data
                css >> i, css >> str, css >> str, css >> x, css >> y, css >> z;

                // set the data
                G(i - 1, 0) = x; G(i - 1, 1) = y; G(i - 1, 2) = z;
            } return G;
        } else if (line == "The final MP2 gradient") {
            // loop over gradient lines
            while (std::getline(lss, line) && !line.empty()) {
                // create the column stringstream
                std::stringstream css; css << line;

                // exctract data
                css >> str, css >> x, css >> y, css >> z;

                // get the row index
                int i = std::stoi(str.substr(0, str.size() - 1));

                // set the data
                G(i, 0) = x; G(i, 1) = y; G(i, 2) = z;
            } return G;
        }
    }

    // return the gradient
    return G;
}
