#pragma once

#include "argparse.hpp"
#include "constant.h"

#include <unordered_map>
#include <filesystem>

class Parser {
public:
    // constructors
    Parser(const argparse::ArgumentParser& program) : program(program) {} Parser(int argc, char** argv);
    Parser(const std::string& name) : program(name, "0.1", argparse::default_arguments::none) {}

    // methods
    bool used(const std::string& name) const {return program.is_subcommand_used(name);}
    bool has(const std::string& name) const {return program.is_used(name);}
    Parser& at(const std::string& name) {return parsers.at(name);}
    template <typename T> T get(const std::string& name) const;

public:
    std::unordered_map<std::string, Parser> parsers;
    argparse::ArgumentParser program;
};

template <typename T = std::string>
T Parser::get(const std::string& name) const {
    return program.get<T>(name);
}
