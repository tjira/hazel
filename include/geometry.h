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
        glm::mat4 getModel() const {
            return translate * rotate * scale;
        }
        glm::vec3 getPosition() const {
            return glm::vec3(translate[3]);
        }
        glm::mat4 translate, rotate, scale;
        std::string name;
    };

public:

    // Constructors
    Geometry() {};

    // Statc constructors
    static Geometry Load(std::stringstream& file);

    // Getters
    glm::vec3 getCenter() const;
    size_t size() const;

    // State functions
    void moveBy(const glm::vec3& vector);
    void render(const Shader& shader, const glm::mat4& transform = glm::mat4(1.0f)) const;
    void rebind(float factor);

    // Public static variables
    inline static std::unordered_map<std::string, Mesh> meshes;

private:
    std::vector<Object> objects;
};
