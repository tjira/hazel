#pragma once

#include "glfwpointer.h"
#include "mesh.h"
#include "ptable.h"
#include <glm/gtc/matrix_transform.hpp>
#include <fstream>
#include <sstream>
#include <unordered_map>

class Geometry {

    // Private subclasses
    struct Object {
        glm::vec3 getPosition() const {
            return glm::vec3(model[3]);
        }
        std::string name;
        glm::mat4 model;
    };

public:

    // Constructors
    Geometry() {};

    // Statc constructors
    static Geometry Load(std::stringstream& file);

    // Getters
    size_t size() const;

    // State functions
    void render(const Shader& shader, const glm::mat4& transform = glm::mat4(1.0f)) const;
    void rebind(float factor);

    // Public static variables
    inline static std::unordered_map<std::string, Mesh> meshes;

private:
    std::vector<Object> objects;
};
