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
        input["bagel"].push_back("{\"title\":\"mp2\",\"aux_basis\":\"cc-pvdz-ri\"}"_json);
    } else if (Utility::StringContains(opt.method, "fci")) {
        input["bagel"].push_back("{\"title\":\"hf\"}"_json);
        input["bagel"].push_back("{\"title\":\"fci\"}"_json);
    } else if (Utility::StringContains(opt.method, "casscf") || Utility::StringContains(opt.method, "caspt2")) {
        input["bagel"].push_back("{\"title\":\"casscf\"}"_json);
        std::vector<std::string> casopt; size_t last = 0, next = 0;
        while ((next = opt.method.find("/", last)) != std::string::npos) {
            casopt.push_back(opt.method.substr(last, next - last)); last = next + 1;
        } casopt.push_back(opt.method.substr(last));
        input["bagel"][1]["nstate"] = std::stoi(casopt.at(3)), input["bagel"][1]["nact"] = std::stoi(casopt.at(2));
        input["bagel"][1]["nclosed"] = system.electrons / 2 - std::stoi(casopt.at(1)) / 2;
        input["bagel"][1]["maxiter"] = 1000, input["bagel"][1]["maxiter_micro"] = 1000;
        if (Utility::StringContains(opt.method, "caspt2")) {
            input["bagel"].push_back("{\"title\":\"smith\"}"_json);
            input["bagel"][2]["method"] = "caspt2";
            input["bagel"][2]["shift"] = 0.2;
        }
    } else {
        throw std::runtime_error("METHOD NOT IMPLEMENTED IN HAZEL-BAGEL INTERFACE");
    }

    // write the coordinates
    for (int i = 0; i < system.coords.rows(); i++) {
        input["bagel"][0]["geometry"].push_back({{"atom", an2sm.at(system.atoms.at(i).atomic_number)}, {"xyz", {system.coords(i, 0), system.coords(i, 1), system.coords(i, 2)}}});
    }
}

void Bagel::enableGradient(const std::vector<int>& targets) {
    // create gradient methods
    nlohmann::json gradient = {{"title", "force"}, {"method", {}}, {"export", true}};

    // add gradient targets
    if (input["bagel"][input["bagel"].size() - 1]["title"] == "casscf" || input["bagel"][input["bagel"].size() - 1]["title"] == "smith") {
        // rename the block title
        gradient["title"] = "forces";

        // push the targets
        for (int state : targets) {
            gradient["grads"].push_back({{"title", "force"}, {"target", state}});
        }
    }

    // fill the gradient block
    if (input["bagel"][input["bagel"].size() - 1]["title"] == "smith") {
        gradient["method"] = {{
            {"title", "caspt2"},
            {"nclosed", input["bagel"][1]["nclosed"]},
            {"nstate", input["bagel"][1]["nstate"]},
            {"nact", input["bagel"][1]["nact"]},
            {"smith", {
                {"method", "caspt2"},
                {"shift", 0.2}
            }}
        }};
    } else {
        for (size_t i = 1; i < input["bagel"].size(); i++) {
            gradient["method"].push_back(input["bagel"][i]);
        }
    }

    // save the targets
    this->targets = targets;

    // modify input
    input["bagel"].push_back(gradient);
}

void Bagel::enableHessian(double) {
    // create hessian methods
    nlohmann::json hessian = {{"title", "hessian"}, {"method", {}}, {"export", true}};

    // fill the gradient block
    if (input["bagel"][input["bagel"].size() - 1]["title"] == "smith") {
        hessian["method"] = {{
            {"title", "caspt2"},
            {"nclosed", input["bagel"][1]["nclosed"]},
            {"nstate", input["bagel"][1]["nstate"]},
            {"nact", input["bagel"][1]["nact"]},
            {"smith", {
                {"method", "caspt2"},
                {"shift", 0.2}
            }}
        }};
    } else {
        for (size_t i = 1; i < input["bagel"].size(); i++) {
            hessian["method"].push_back(input["bagel"][i]);
        }
    }

    // modify input
    input["bagel"].push_back(hessian);
}

Bagel::Results Bagel::run() const {
    // create directory and open the file
    std::filesystem::create_directory(directory);
    std::ofstream ifile(directory + "/bagel.json");

    // write the input to the opened file
    ifile << input.dump(2) << std::endl; ifile.close();

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
    std::vector<double> excs; std::smatch match; int nstate = 0;
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
    return (excs.size() ? Eigen::Map<Vector>(excs.data() + excs.size() - nstate, nstate) : Vector());
}

double Bagel::extractEnergy(const std::string& output) const {
    // define iterator, smatch and the placeholder for the energy
    std::string::const_iterator sstart(output.begin());
    std::smatch match; double energy = 0;

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
    std::string::const_iterator sstart(output.begin());
    std::vector<double> f; std::smatch match;

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

std::vector<Matrix> Bagel::extractGradient(const std::string&) const {
    // create the gradient vector
    std::vector<Matrix> Gs;

    for (int i : targets) {
        // create the gradient matrix
        Matrix G(system.coords.rows(), 3); std::string line;
        std::ifstream fstream(directory + "/FORCE_" + std::to_string(i) + ".out");

        // loop over lines of gradient
        for (int i = 0; std::getline(fstream, line); i++) if (i && line.size()) {
            // creale column stringstream
            std::stringstream css(line); css >> i;

            // assign gradient values
            css >> G(i, 0), css >> G(i, 1), css >> G(i, 2);
        }

        // append the gradient
        Gs.push_back(G);
    }

    // return the gradient
    return Gs;
}
