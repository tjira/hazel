#include "../include/movie.h"

Movie Movie::LoadTrajectory(const std::string& filename) {
    Movie movie;

    // Open a file, create buffer and clear previous molecule.
    std::ifstream file(filename); std::string line;
    std::getline(file, line);

    // Get the atom count (length) and file length (block).
    int length = std::stoi(line), block = 1;
    while (std::getline(file, line)) block++;

    // Create the vector of geometris and reset the file reader.
    movie.scenes.resize(block / (length + 2));
    file.clear(), file.seekg(0);

    // Read the individual geometries.
    for (int i = 0; i < block / (length + 2); i++) {
        std::stringstream ss;

        // For lines in one geometry.
        for (int j = 0; j < length + 2; j++) {
            std::getline(file, line); ss << (j ? "\n" : "" ) << line;
        }

        // Load a geometry to a molecule class.
        movie.scenes.at(i) = Scene::LoadMolecule(ss);
    }

    // Set the initialization timestamp (for FPS manipulation).
    movie.timestamp = std::chrono::high_resolution_clock().now();

    return movie;
}

void Movie::render(const Shader& shader) {
    if (scenes.size()) {
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock().now() - timestamp).count();
        if (elapsed > 16) {
            frame = (frame + (int)(elapsed / 16)) % (int)scenes.size();
            timestamp = std::chrono::high_resolution_clock().now();
        }
        scenes.at(frame).render(shader);
    }
}
