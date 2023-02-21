#pragma once

#include "moleculegraphic.h"

class TrajectoryGraphic {
public:
    
    // Constructors
    TrajectoryGraphic() {}

    // Static constructors
    static TrajectoryGraphic Load(const std::string& movie);

    // Getters
    std::vector<MoleculeGraphic>& getGeoms() { return geoms; }
    bool& getPause() { return paused; }
    int& getFrame() { return frame; }
    int size() const { return geoms.size(); }

    // State functions
    void render(const Shader& shader);

private:
    std::chrono::high_resolution_clock::time_point timestamp;
    std::vector<MoleculeGraphic> geoms;
    bool paused = false;
    int frame = 0;
};
