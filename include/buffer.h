#pragma once

#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <vector>

struct Vertex {
    glm::vec3 position, normal, color = glm::vec3(1);
};

class Buffer {
public:
    Buffer(const Buffer& buffer) : data(buffer.getData()) { generate(); };
    Buffer(const std::vector<Vertex>& data) : data(data) { generate(); };
    Buffer() { generate(); };
    ~Buffer() { glDeleteVertexArrays(1, &vao), glDeleteBuffers(1, &vbo); };
    Buffer& operator=(const Buffer& buffer);
    void bind() const { glBindVertexArray(vao); };
    std::vector<Vertex> getData() const { return data; }
    size_t getSize() const { return data.size(); };

private:
    void generate();
    unsigned int vao, vbo;
    std::vector<Vertex> data;
};
