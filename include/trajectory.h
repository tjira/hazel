#pragma once

#include "geometry.h"

class Trajectory {
public:
    
    // Constructors
    Trajectory() {}

    // Static constructors
    static Trajectory Load(const std::string& movie);

    // Getters
    std::vector<Geometry>& getGeoms() { return geoms; }
    bool& getPause() { return paused; }
    int& getFrame() { return frame; }
    float& getWait() { return wait; }
    int size() const { return geoms.size(); }

    // State functions
    void moveBy(const glm::vec3& vector);
    void render(const Shader& shader);

private:
    std::chrono::high_resolution_clock::time_point timestamp;
    std::vector<Geometry> geoms;
    bool paused = false;
    float wait = 16.0;
    int frame = 0;
};
