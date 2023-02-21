#include "../include/scene.h"

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
        scene.objects.push_back({ atom, glm::scale(glm::translate(glm::mat4(1.0f), glm::vec3(x, y, z)), glm::vec3(0.007f * ptable.at(atom).radius)) });
     }

    // Add bonds
    scene.rebindMolecule(BINDINGFACTOR);

    // Return scene
    return scene;
}

void Scene::rebindMolecule(float factor) {
    std::vector<Object> objects;
    for (Object obj : this->objects) {
        if (obj.name != "bond") objects.push_back(obj);
    }

    int length = objects.size();

    for (size_t i = 0; i < length; i++) {
        for (size_t j = i + 1; j < length; j++) {
            float distance = glm::length(objects.at(j).getPosition() - objects.at(i).getPosition());
            if (objects.at(i).name == "El" || objects.at(j).name == "El") continue;
            if (distance < factor * (ptable.at(objects.at(i).name).covalent + ptable.at(objects.at(j).name).covalent)) {
                glm::vec3 position = (objects.at(i).getPosition() + objects.at(j).getPosition()) / 2.0f;
                glm::vec3 vector = objects.at(j).getPosition() - objects.at(i).getPosition();
                glm::vec3 cross = glm::cross(glm::vec3(0, 1, 0), vector);
                float angle = atan2f(glm::length(cross), glm::dot(glm::vec3(0, 1, 0), vector));
                glm::mat4 model = glm::scale(glm::rotate(glm::translate(glm::mat4(1.0f), position), angle, glm::normalize(cross)), { 0.09f, glm::length(vector) / 2.0f, 0.09f });
                objects.push_back({ "bond", model });
            }
        }
    }

    this->objects = objects;
};

void Scene::render(const Shader& shader, const glm::mat4& transform) const {
    for (size_t i = 0; i < objects.size(); i++) {
        if (objects.at(i).name == "bond") meshes.at("bond").render(shader, transform * objects.at(i).model);
        else meshes.at(objects.at(i).name).render(shader, transform * objects.at(i).model);
    }
}

size_t Scene::size() const {
    return objects.size();
}
