#pragma once

#include "scene.h"

class Movie {
public:
    
    // Constructors
    Movie() {}

    // Static constructors
    static Movie LoadTrajectory(const std::string& movie);

    // State functions
    void render(const Shader& shader);

private:
    std::chrono::high_resolution_clock::time_point timestamp;
    std::vector<Scene> scenes;
    int frame = 0;
};
