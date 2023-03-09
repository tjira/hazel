#pragma once

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace Defaults {
    static json hfopt = R"({
        "name" : "HF",
        "maxiter" : 100,
        "thresh" : 1e-8,
        "diis" : {
            "start" : 3,
            "keep" : 5,
            "damp" : 0
        },
        "engrad" : {
            "numerical" : true,
            "increment" : 0.005,
            "nthread" : 4
        },
        "mulliken" : false,
        "print" : {
            "kinetic" : false,
            "oneelec" : false,
            "overlap" : false,
            "orben" : false,
            "density" : false
        }
    })"_json;
    static json mdopt = R"({
        "timestep" : 0.001,
        "steps" : 10000,
        "output" : "trajectory.xyz"
    })"_json;
}
