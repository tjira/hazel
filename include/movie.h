#pragma once

#include "scene.h"

class Movie {
public:
    
    // Constructors
    Movie() {}

    // Static constructors
    static Movie LoadTrajectory(const std::string& movie);

    // Getters
    bool& getPause() { return paused; }
    int& getFrame() { return frame; }
    int size() const { return scenes.size(); }

    // State functions
    void render(const Shader& shader);

private:
    std::chrono::high_resolution_clock::time_point timestamp;
    std::vector<Scene> scenes;
    bool paused = false;
    int frame = 0;
};
