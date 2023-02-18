#include "../include/scene.h"

Scene Scene::LoadMolecule(const std::string& filename) {
    std::ifstream file(filename);
    std::stringstream ss;
    ss << file.rdbuf();
    return LoadMolecule(ss);
}

Scene Scene::LoadMolecule(std::stringstream& file) {
    // Open file and declare variables
    Scene scene; int length;
    std::string line, name;

    // Extract length and name.
    std::getline(file, line);
    std::getline(file, name);
    length = std::stoi(line);

    // Get first line of coordinates
    std::getline(file, line);

    // Add atom for each line.
    for (int i = 0; i < length; std::getline(file, line), i++) {
        std::string atom; float x, y, z;
        std::stringstream iss(line);
        iss >> atom >> x >> y >> z;
        scene.objects.push_back(Mesh::Icosphere(1, 0, atom));
        scene.at(i).setModel(glm::scale(glm::translate(glm::mat4(1.0f), glm::vec3(x, y, z)), glm::vec3(0.007f * ptable.at(atom).radius))); 
        scene.at(i).setColor(ptable.at(atom).color);
     }

    // Add bonds
    for (size_t i = 0; i < length; i++) {
        for (size_t j = i + 1; j < length; j++) {
            float distance = glm::length(scene.objects.at(j).getPosition() - scene.objects.at(i).getPosition());
            if (distance < 0.013f * (ptable.at(scene.objects.at(i).getName()).covalent + ptable.at(scene.objects.at(j).getName()).covalent)) {
                glm::vec3 position = (scene.objects.at(i).getPosition() + scene.objects.at(j).getPosition()) / 2.0f;
                glm::vec3 vector = scene.objects.at(j).getPosition() - scene.objects.at(i).getPosition();
                glm::vec3 cross = glm::cross(glm::vec3(0, 1, 0), vector);
                float angle = atan2f(glm::length(cross), glm::dot(glm::vec3(0, 1, 0), vector));
                glm::mat4 model = glm::scale(glm::rotate(glm::translate(glm::mat4(1.0f), position), angle, glm::normalize(cross)), { 0.09f, glm::length(vector) / 2.0f, 0.09f });
                scene.objects.push_back(Mesh::Cylinder(8, 0, "bond"));
                scene.at(scene.size() - 1).setModel(model);
            }
        }
    }
    return scene;
}

Mesh& Scene::at(int i) {
    return objects.at(i);
}

void Scene::render(const Shader& shader, const glm::mat4& transform) const {
    for (size_t i = 0; i < size(); i++) {
        objects.at(i).render(shader, transform);
    }
}

size_t Scene::size() const {
    return objects.size();
}
