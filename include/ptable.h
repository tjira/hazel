#pragma once

#include <unordered_map>

struct Atom {
    double mass;
};

namespace {
    inline static std::unordered_map<int, Atom> ptable;
}
