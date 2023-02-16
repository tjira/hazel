#pragma once

#include "buffer.h"
#include "shader.h"
#include <filesystem>
#include <functional>
#include <iostream>
#include <stdexcept>

class Mesh {
public:
    Mesh() {}; Mesh(std::vector<Vertex3D> data, const std::string& name = "mesh") : name(name), model(1.0f), buffer(data) {};
    static Mesh Cylinder(int sectors, bool smooth, const std::string& name = "cylinder");
    static Mesh Icosphere(int subdivisions, bool smooth, const std::string& name = "icosphere");
    void render(const Shader& shader, const glm::mat4& transform = glm::mat4(1.0f)) const;

private:
    std::string name;
    glm::mat4 model;
    Buffer buffer;
};
