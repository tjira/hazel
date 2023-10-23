#pragma once

#include "argparse.hpp"
#include "constant.h"

#include <unordered_map>
#include <filesystem>

class Parser {
public:
    // constructors
    Parser(const std::string& name) : program(name, "0.1", argparse::default_arguments::none), name(name) {}
    Parser(const argparse::ArgumentParser& program) : program(program) {} Parser(int argc, char** argv);

    // methods
    bool used(const std::string& name) const {return program.is_subcommand_used(name);}
    bool has(const std::string& name) const {return program.is_used(name);}
    Parser& at(const std::string& name) {return *parsers.at(name);}
    template <typename T> T get(const std::string& name) const;
    void help() const;

public:
    std::unordered_map<std::string, std::shared_ptr<Parser>> parsers;
    argparse::ArgumentParser program; std::string name;
};

template <typename T>
T Parser::get(const std::string& name) const {
    return program.get<T>(name);
}
