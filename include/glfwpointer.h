#pragma once

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <string>

#define WIDTH 1024
#define HEIGHT 576

struct GLFWPointer {
    std::string title = "Hazel Viewer"; glm::vec2 mouse; GLFWwindow* window;
    int width = WIDTH, height = HEIGHT, samples = 16, major = 4, minor = 2;
    struct Camera {
        glm::mat4 view, proj;
    } camera{};
    struct Flags {
        bool fullscreen = false, info = false, renderOptions = false;
    } flags{};
    struct Light {
        float ambient = 0.4f, diffuse = 0.2f, specular = 0.4f, shininess = 4.0f;
        glm::vec3 position = { 1.0f, 1.0f, 1.0f };
    } light{};
};
