#pragma once

#include "mesh.h"
#include "ptable.h"
#include <glm/gtc/matrix_transform.hpp>
#include <fstream>
#include <sstream>

class Scene {
public:
    Scene() {};
    static Scene LoadMolecule(std::stringstream& file);
    Mesh& at(int i);
    size_t size() const;
    void render(const Shader& shader, const glm::mat4& transform = glm::mat4(1.0f)) const;

private:
    std::vector<Mesh> objects;
};
