#include "orca.h"

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
        input = std::regex_replace(input, std::regex("\n\\*xyz"), "\n%casscf\nnel " + casopt.at(1) + "\nnorb " + casopt.at(2) + "\nmult 1\nnroots " + casopt.at(3) + "\nend\n\n*xyz");
    } else {
        input = std::regex_replace(input, std::regex("METHOD"), Utility::ToUpper(opt.method) + " HCORE");
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

    // write the input to the opened file
    ifile << input; ifile.close();

    // define the buffer for output
    std::array<char, 128> buffer; std::string output;

    // execute the command
    #ifdef _WIN32
    auto pipe = _popen(("cd " + directory + " && orca orca.inp > >(tee orca.out) 2> /dev/null").c_str(), "r");
    #else
    auto pipe = popen(("cd " + directory + " && orca orca.inp > >(tee orca.out) 2> /dev/null").c_str(), "r");
    #endif

    // check for success
    if (!pipe) throw std::runtime_error("ORCA EXECUTION FAILED");

    // read the output
    while (!feof(pipe)) {
        if (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
            output += buffer.data();
        }
    }

    // throw error
    #ifdef _WIN32
    if (_pclose(pipe) != EXIT_SUCCESS) {
        throw std::runtime_error("ORCA TERMINATED UNSUCCESSFULLY");
    }
    #else
    if (pclose(pipe) != EXIT_SUCCESS) {
        throw std::runtime_error("ORCA TERMINATED UNSUCCESSFULLY");
    }
    #endif

    // copy the output
    std::filesystem::copy_file(directory + "/orca.out", "orca.out", std::filesystem::copy_options::overwrite_existing);

    // remove the calculataion directory
    std::filesystem::remove_all(directory);

    // extract the results
    Results results = {extractEnergy(output), extractEnergies(output), extractFrequencies(output), extractGradient(output)};

    // assign the correct ground state energy
    if (results.excs.size()) results.E = results.excs(0);

    // return the results
    return results;
}

Vector Orca::extractEnergies(const std::string& output) const {
    // define the start iterator and match
    std::string::const_iterator sstart(output.begin());
    std::smatch match; std::vector<double> excs;

    // find the number of roots
    if (std::regex_search(output, match, std::regex(".* NROOTS\\=\\ *(\\d*)"))) {
        // assign the number of roots
        int nroots = std::stoi(match[1]);

        // find append the CASSCF roots to the result
        while (std::regex_search(sstart, output.end(), match, std::regex(".* E\\=\\ *([\\-\\.\\d]*) Eh.*"))) {
            excs.push_back(std::stod(match[1])), sstart = match.suffix().first;
        }

        // return the results
        if (excs.size() >= (size_t)nroots) return Eigen::Map<Vector>(excs.data() + excs.size() - nroots, nroots);
    }

    // return nothing
    return Vector();
}

double Orca::extractEnergy(const std::string& output) const {
    // results in the following order: FCI, OTHER, NOT FOUND
    if (std::smatch match; std::regex_search(output, match, std::regex("FINAL SINGLE .* Root 0 \t\\= ([\\-\\.\\d]*)"))) {
        return std::stod(match[1]);
    } else if (std::regex_search(output, match, std::regex("FINAL SINGLE .* ([\\-\\.\\d]*)"))) {
        return std::stod(match[1]);
    } else return 0;
}

Vector Orca::extractFrequencies(const std::string& output) const {
    // define match, iterator and frequency vector
    std::string::const_iterator sstart(output.begin());
    std::smatch match; std::vector<double> f;

    // find append the CASSCF roots to the result
    while (std::regex_search(sstart, output.end(), match, std::regex(".* ([\\-\\.\\d]+) cm.*"))) {
        f.insert(f.begin(), std::stod(match[1])), sstart = match.suffix().first;
    }

    // return the frequencies
    return (f.size() ? Eigen::Map<Vector>(f.data(), f.size()) : Vector());
}

Matrix Orca::extractGradient(const std::string& output) const {
    // define match, iterator and gradient matrix
    std::smatch match; Matrix G(system.coords.rows(), 3);
    std::string::const_iterator sstart(output.begin());

    // find append the CASSCF roots to the result
    for (int i = 0; std::regex_search(sstart, output.end(), match, std::regex(".* \\d+.*:\\ *([\\-\\.\\d]+)\\ +([\\-\\.\\d]+)\\ +([\\-\\.\\d]+)\\n")); i++) {
        G(i, 0) = std::stod(match[1]), G(i, 1) = std::stod(match[2]), G(i, 2) = std::stod(match[3]), sstart = match.suffix().first;
    }

    // return gradient
    return G;
}
