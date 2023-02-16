#pragma once

#include <glm/glm.hpp>
#include <unordered_map>

struct Atom {
    float radius, covalent;
    glm::vec3 color;
    double mass;
};

extern std::unordered_map<std::string, Atom> ptable;
