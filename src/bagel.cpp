#include "bagel.h"

Bagel::Bagel(const System& system, const Options& opt) : opt(opt), system(system), input(nlohmann::json::parse(BAGEL)) {
    // define the directory name and remove the starting new line
    directory = ".bagel." + std::to_string(Timer::Now().time_since_epoch().count());

    // set basis
    input["bagel"][0]["basis"] = system.originbasis;

    // replace placeholders for method
    if (Utility::StringContains(opt.method, "hf")) {
        input["bagel"].push_back("{\"title\":\"hf\"}"_json);
    } else if (Utility::StringContains(opt.method, "mp2")) {
        input["bagel"].push_back("{\"title\":\"hf\"}"_json);
        input["bagel"].push_back("{\"title\":\"mp2\",\"aux_basis\":\"cc-pvdz-ri\"}"_json);
    } else if (Utility::StringContains(opt.method, "fci")) {
        input["bagel"].push_back("{\"title\":\"hf\"}"_json);
        input["bagel"].push_back("{\"title\":\"fci\"}"_json);
        std::vector<std::string> casopt; size_t last = 0, next = 0;
        while ((next = opt.method.find("/", last)) != std::string::npos) {
            casopt.push_back(opt.method.substr(last, next - last)); last = next + 1;
        } casopt.push_back(opt.method.substr(last));
        input["bagel"][2]["nstate"] = std::stoi(casopt.at(1));
    } else if (Utility::StringContains(opt.method, "casscf") || Utility::StringContains(opt.method, "caspt2")) {
        input["bagel"].push_back("{\"title\":\"casscf\"}"_json);
        std::vector<std::string> casopt; size_t last = 0, next = 0;
        while ((next = opt.method.find("/", last)) != std::string::npos) {
            casopt.push_back(opt.method.substr(last, next - last)); last = next + 1;
        } casopt.push_back(opt.method.substr(last));
        input["bagel"][1]["nstate"] = std::stoi(casopt.at(3)), input["bagel"][1]["nact"] = std::stoi(casopt.at(2));
        input["bagel"][1]["nclosed"] = system.electrons / 2 - std::stoi(casopt.at(1)) / 2;
        if (Utility::StringContains(opt.method, "caspt2")) {
            input["bagel"].push_back("{\"title\":\"smith\"}"_json);
            input["bagel"][2]["method"] = "caspt2";
            input["bagel"][2]["ms"] = true;
            input["bagel"][2]["xms"] = true;
            input["bagel"][2]["sssr"] = true;
            input["bagel"][2]["shift"] = 0.2;
            input["bagel"][2]["nstate"] = std::stoi(casopt.at(3)), input["bagel"][2]["nact"] = std::stoi(casopt.at(2));
            input["bagel"][2]["nclosed"] = system.electrons / 2 - std::stoi(casopt.at(1)) / 2;
        }
    } else {
        throw std::runtime_error("METHOD NOT IMPLEMENTED IN HAZEL-BAGEL INTERFACE");
    }

    // write the coordinates
    for (int i = 0; i < system.coords.rows(); i++) {
        input["bagel"][0]["geometry"].push_back({{"atom", an2sm.at(system.atoms.at(i).atomic_number)}, {"xyz", {system.coords(i, 0), system.coords(i, 1), system.coords(i, 2)}}});
    }
}

void Bagel::enableGradient(double) {
    input["bagel"][input["bagel"].size() - 1] = {{"title", "force"}, {"method", {input["bagel"][input["bagel"].size() - 1]}}};
    input["bagel"][input["bagel"].size() - 1]["export"] = true;
    if (input["bagel"][input["bagel"].size() - 1]["method"][0]["title"] == "smith") {
        // change title
        input["bagel"][input["bagel"].size() - 1]["method"][0]["title"] = "caspt2";

        // erase
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("method");
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("shift");
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("sssr");
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("xms");
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("ms");

        // smith block
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["method"] = "caspt2";
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["sssr"] = true;
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["shift"] = 0.2;
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["xms"] = true;
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["ms"] = true;
    } else if (input["bagel"][input["bagel"].size() - 1]["method"][0]["title"] == "fci") {
        throw std::runtime_error("COULD NOT CALCULATE GRADIENT FOR FCI");
    } else if (input["bagel"][input["bagel"].size() - 1]["method"][0]["title"] == "hessian") {
        throw std::runtime_error("CALCULATION OF GRADIENT AND HESSIAN NOT POSSIBLE");
    }
}

void Bagel::enableHessian(double) {
    input["bagel"].push_back({{"title", "hessian"}, {"method", {input["bagel"][input["bagel"].size() - 1]}}});
    if (input["bagel"][input["bagel"].size() - 1]["method"][0]["title"] == "smith") {
        // change title
        input["bagel"][input["bagel"].size() - 1]["method"][0]["title"] = "caspt2";

        // erase
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("method");
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("shift");
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("sssr");
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("xms");
        input["bagel"][input["bagel"].size() - 1]["method"][0].erase("ms");

        // smith block
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["method"] = "caspt2";
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["sssr"] = true;
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["shift"] = 0.2;
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["xms"] = true;
        input["bagel"][input["bagel"].size() - 1]["method"][0]["smith"]["ms"] = true;
    } else if (input["bagel"][input["bagel"].size() - 1]["method"][0]["title"] == "fci") {
        throw std::runtime_error("COULD NOT CALCULATE HESSIAN FOR FCI");
    } else if (input["bagel"][input["bagel"].size() - 1]["method"][0]["title"] == "force") {
        throw std::runtime_error("CALCULATION OF GRADIENT AND HESSIAN NOT POSSIBLE");
    }
}

Bagel::Results Bagel::run() const {
    // create directory and open the file
    std::filesystem::create_directory(directory);
    std::ofstream ifile(directory + "/bagel.json");

    // write the input to the opened file
    ifile << input.dump(2); ifile.close();

    // define the buffer for output
    std::array<char, 128> buffer; std::string output;

    // execute the command
    #ifdef _WIN32
    auto pipe = _popen(("cd " + directory + " && BAGEL bagel.json > >(tee bagel.out) 2> /dev/null").c_str(), "r");
    #else
    auto pipe = popen(("cd " + directory + " && BAGEL bagel.json > >(tee bagel.out) 2> /dev/null").c_str(), "r");
    #endif

    // check for success
    if (!pipe) throw std::runtime_error("BAGEL EXECUTION FAILED");

    // read the output
    while (!feof(pipe)) {
        if (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
            output += buffer.data();
        }
    }

    #ifdef _WIN32
    if (_pclose(pipe) != EXIT_SUCCESS) {
        throw std::runtime_error("BAGEL TERMINATED UNSUCCESSFULLY");
    }
    #else
    if (pclose(pipe) != EXIT_SUCCESS) {
        throw std::runtime_error("BAGEL TERMINATED UNSUCCESSFULLY");
    }
    #endif

    // copy the output
    std::filesystem::copy_file(directory + "/bagel.out", "bagel.out", std::filesystem::copy_options::overwrite_existing);
    
    // return the results
    Results results = {extractEnergy(output), extractEnergies(output), extractFrequencies(output), extractGradient(output)};

    // remove the calculataion directory
    std::filesystem::remove_all(directory);

    // assign the correct ground state
    if (results.excs.size()) {
        results.E = results.excs(0);
    }

    // return the results
    return results;
}

Vector Bagel::extractEnergies(const std::string& output) const {
    // define iterator, smatch and number of states
    std::vector<double> excs; std::smatch match; int nstate;
    std::string::const_iterator sstart(output.begin());

    // find and append the FCI or CASSCF roots to the result
    while (std::regex_search(sstart, output.end(), match, std::regex(".*(\\d+) \\*\\ *([\\-\\.\\d]+).*"))) {
        nstate = std::stoi(match[1]) + 1, excs.push_back(std::stod(match[2])), sstart = match.suffix().first;
    }

    // find and append the CAPT2 roots to the result
    while (std::regex_search(sstart, output.end(), match, std::regex(".*MS-CASPT2 energy : state\\ *\\d+\\ *(.*)\n"))) {
        excs.push_back(std::stod(match[1])), sstart = match.suffix().first;
    }

    // return the results
    if (excs.size()) return Eigen::Map<Vector>(excs.data() + excs.size() - nstate, nstate);
    else return Vector();
}

double Bagel::extractEnergy(const std::string& output) const {
    // define iterator, smatch and the placeholder for the energy
    std::string::const_iterator sstart(output.begin());
    std::smatch match; double energy;

    // find append the CASSCF roots to the result
    while (std::regex_search(sstart, output.end(), match, std::regex(".*Fock build.*\n\\ *\\d+\\ *([\\-\\.\\d]+).*"))) {
        energy = std::stod(match[1]), sstart = match.suffix().first;
    }

    if (std::regex_search(sstart, output.end(), match, std::regex("\\ *MP2 total energy:\\ *([\\-\\.\\d]+)"))) {
        energy = std::stod(match[1]), sstart = match.suffix().first;
    }

    // return the results
    return energy;
}

Vector Bagel::extractFrequencies(const std::string& output) const {
    // define iterator, smatch and number of states
    std::vector<double> f; std::smatch match; int nstate;
    std::string::const_iterator sstart(output.begin());

    // find and append the frequencies
    while (std::regex_search(sstart, output.end(), match, std::regex(".*Freq \\(cm\\-1\\)\\ *([\\-\\.\\d]*)\\ *([\\-\\.\\d]*)\\ *([\\-\\.\\d]*)\\ *([\\-\\.\\d]*)\\ *([\\-\\.\\d]*)\\ *([\\-\\.\\d]*)\n"))) {
        for (int i = 1; i < 7 && !match[i].str().empty(); i++) {
            f.insert(f.begin(), std::stod(match[i]));
        }
        sstart = match.suffix().first;
    }

    // return the frequencies
    return (f.size() ? Eigen::Map<Vector>(f.data(), f.size()) : Vector());
}

Matrix Bagel::extractGradient(const std::string& output) const {
    // create the gradient matrix
    Matrix G(system.coords.rows(), 3); std::string line;
    std::ifstream fstream(directory + "/FORCE_0.out");

    // loop over lines of gradient
    for (int i = 0; std::getline(fstream, line); i++) if (i && line.size()) {
        // creale column stringstream
        std::stringstream css(line); css >> i;

        // assign gradient values
        css >> G(i, 0), css >> G(i, 1), css >> G(i, 2);
    }

    // return the gradient
    return G;
}
