#include "bagel.h"

Bagel::Bagel(const System& system, const Options& opt) : opt(opt), system(system), input(nlohmann::json::parse(BAGEL)) {
    // define the directory name and remove the starting new line
    directory = ".bagel." + std::to_string(Timer::Now().time_since_epoch().count());

    // set basis
    input["bagel"][0]["basis"] = system.originbasis;

    // replace placeholders for method
    if (Utility::StringContains(opt.method, "hf")) {
        input["bagel"].push_back("{\"title\":\"hf\"}"_json);
    } else if (Utility::StringContains(opt.method, "casscf")) {
        input["bagel"].push_back("{\"title\":\"casscf\"}"_json);
        std::vector<std::string> casopt; size_t last = 0, next = 0;
        while ((next = opt.method.find("/", last)) != std::string::npos) {
            casopt.push_back(opt.method.substr(last, next - last)); last = next + 1;
        } casopt.push_back(opt.method.substr(last));
        input["bagel"][1]["nstate"] = std::stoi(casopt.at(3)), input["bagel"][1]["nact"] = std::stoi(casopt.at(2));
        input["bagel"][1]["nclosed"] = system.electrons / 2 - std::stoi(casopt.at(1)) / 2;
    } else {
        throw std::runtime_error("METHOD NOT IMPLEMENTED IN HAZEL-BAGEL INTERFACE");
    }

    // write the coordinates
    for (int i = 0; i < system.coords.rows(); i++) {
        input["bagel"][0]["geometry"].push_back({{"atom", an2sm.at(system.atoms.at(i).atomic_number)}, {"xyz", {system.coords(i, 0), system.coords(i, 1), system.coords(i, 2)}}});
    }
}

void Bagel::enableGradient(double) {
    input["bagel"][1] = {{"title", "force"}, {"method", {input["bagel"][1]}}};
}

void Bagel::enableHessian(double) {
    input["bagel"].push_back({{"title", "hessian"}, {"method", {input["bagel"][1]}}});
}

Bagel::Results Bagel::run() const {
    // create directory and open the file
    std::filesystem::create_directory(directory);
    std::ofstream ifile(directory + "/bagel.json");

    // write the input to the opened file
    ifile << input.dump(2); ifile.close();

    // define the buffer for output and result string with struct
    std::array<char, 128> buffer; std::string output; Results results;

    // execute the command
    #ifdef _WIN32
    std::unique_ptr<FILE, decltype(&pclose)> pipe(_popen(("cd " + directory).c_str(), "r"), _pclose);
    #else
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(("cd " + directory + " && BAGEL bagel.json > >(tee bagel.out) 2> /dev/null").c_str(), "r"), pclose);
    #endif

    // check for success
    if (!pipe) throw std::runtime_error("BAGEL EXECUTION FAILED");

    // read the output
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        if (Utility::StringContains(std::string(buffer.data()), "Error:")) {
            throw std::runtime_error(std::string(buffer.data()));
        } else if (Utility::StringContains(std::string(buffer.data()), "ERROR:")) {
            throw std::runtime_error(std::string(buffer.data()));
        } else if (Utility::StringContains(std::string(buffer.data()), "INPUT ERROR")) {
            throw std::runtime_error(std::string(buffer.data()));
        }
        output += std::string(buffer.data());
    }

    // copy the output
    std::filesystem::copy_file(directory + "/bagel.out", "bagel.out", std::filesystem::copy_options::overwrite_existing);

    // remove the calculataion directory
    std::filesystem::remove_all(directory);
    
    // return the results
    return {extractEnergy(output), extractEnergies(output), extractFrequencies(output), extractGradient(output)};
}

Vector Bagel::extractEnergies(const std::string& output) const {
    // create the string stream and line buffer
    std::stringstream lss; lss << output; std::string line;

    // create the energy vector and some variables
    std::vector<double> excs; std::string str; double energy; int states;

    // loop over lines
    while (std::getline(lss, line)) {
        if (Utility::StringContains(line, "nstate")) {std::stringstream css; css << line; css >> str, css >> str, css >> str, css >> states;}
        if (Utility::StringContains(line, "FCI iteration")) {
            // loop over excitation lines
            while (std::getline(lss, line) && !Utility::StringContains(line, "vector")) {
                if (Utility::StringContains(line, "*")) {
                    // create the column stringstream
                    std::stringstream css; css << line;

                    // exctract the correct column energy
                    for (int i = 0; i < 3; i++) {css >> str;} css >> energy;

                    // append the energy
                    excs.push_back(energy);
                }
            } return Eigen::Map<Vector>(excs.data() + excs.size() - states, states);
        }
    }

    // return the results
    return Vector();
}

double Bagel::extractEnergy(const std::string& output) const {
    // create the string stream, line buffer and energy placeholder
    std::stringstream lss; lss << output; std::string line; double E = 0;

    // loop over lines
    while (std::getline(lss, line)) {
        if (Utility::StringContains(line, "Fock build")) {
            // skil line and create the column stringstream
            std::getline(lss, line); std::stringstream css; css << line;

            // exctract the correct column energy
            std::string cell; for (int i = 0; i < 2; i++) css >> cell;

            // assign the final energy
            E = std::stod(cell);
        }
    }

    // return the results
    return E;
}

Vector Bagel::extractFrequencies(const std::string& output) const {
    // create the string stream and line buffer
    std::stringstream lss; lss << output; double freq;
    std::string line, str; std::vector<double> f;

    // loop over lines
    while (std::getline(lss, line)) {
        if (Utility::StringContains(line, "Freq")) {
            std::stringstream css; css << line;

            // skip columns
            for (int i = 0; i < 2; i++) css >> str;

            // for every frequency
            while (css >> freq) {
                if (freq) f.push_back(freq);
            }
        }
    }

    // return the frequencies
    std::reverse(f.begin(), f.end()); return Eigen::Map<Vector>(f.data(), f.size());
}

Matrix Bagel::extractGradient(const std::string& output) const {
    // create the string stream and line buffer
    std::stringstream lss; lss << output; std::string line;

    // create the gradient matrix and some variables
    Matrix G = Matrix::Zero(system.coords.rows(), 3);
    double val; std::string str;

    // loop over lines
    while (std::getline(lss, line)) {
        if (Utility::StringContains(line, "Nuclear energy gradient")) {
            // skip lines
            for (int i = 0; i < 1; i++) std::getline(lss, line);

            // loop over gradient lines
            for (size_t i = 0; i < system.atoms.size(); i++) {
                std::getline(lss, line);

                // loop over gradient coordinates
                for (int j = 0; j < 3; j++) {
                    std::getline(lss, line); std::stringstream css; css << line;
                    css >> str, css >> val; G(i, j) = val;
                }

                // set the data
            } return G;
        }
    }

    // return the gradient
    return G;
}
