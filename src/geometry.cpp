#include "../include/geometry.h"

/*
Read the geometry from an .xyz file.
*/
Geometry Geometry::Load(std::stringstream& file) {
    // Open file and declare variables
    Geometry molecule; int length;
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
        glm::mat4 scale = glm::scale(glm::mat4(1), glm::vec3(0.007f * ptable.at(atom).radius));
        glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(x, y, z));
        molecule.objects.push_back({ translate, glm::mat4(1.0f), scale, atom });
     }

    // Add bonds
    molecule.rebind(BINDINGFACTOR);

    // Return molecule
    return molecule;
}

/*
Returns the geometric center of the molecule.
*/
glm::vec3 Geometry::getCenter() const {
    glm::vec3 center(0); float size = 0;
    for (const Object& object : objects) {
        if (object.name != "bond" && object.name != "El") {
            center += object.getPosition(), size += 1;
        }
    }
    return center / size;
}

/*
Move the molecule by some vector.
*/
void Geometry::moveBy(const glm::vec3& vector) {
    for (Object& object : objects) {
        object.translate = glm::translate(object.translate, vector);
    }
}

/*
Create bonds for atoms based on the binding factor.
*/
void Geometry::rebind(float factor) {
    std::vector<Object> objects;
    for (Object obj : this->objects) {
        if (obj.name != "bond") objects.push_back(obj);
    }

    size_t length = objects.size();

    for (size_t i = 0; i < length; i++) {
        for (size_t j = i + 1; j < length; j++) {
            float distance = glm::length(objects.at(j).getPosition() - objects.at(i).getPosition());
            if (objects.at(i).name == "El" || objects.at(j).name == "El") continue;
            if (distance < factor * (ptable.at(objects.at(i).name).covalent + ptable.at(objects.at(j).name).covalent)) {
                glm::vec3 position = (objects.at(i).getPosition() + objects.at(j).getPosition()) / 2.0f;
                glm::vec3 vector = objects.at(j).getPosition() - objects.at(i).getPosition();
                glm::vec3 cross = glm::cross(glm::vec3(0, 1, 0), vector);
                float angle = atan2f(glm::length(cross), glm::dot(glm::vec3(0, 1, 0), vector));
                glm::mat4 scale = glm::scale(glm::mat4(1), { 0.09f, glm::length(vector) / 2.0f, 0.09f });
                glm::mat4 rotate = glm::rotate(glm::mat4(1), angle, glm::normalize(cross));
                glm::mat4 translate = glm::translate(glm::mat4(1.0f), position);
                objects.push_back({ translate, rotate, scale, "bond" });
            }
        }
    }

    this->objects = objects;
};

/*
Render the geometry
*/
void Geometry::render(const Shader& shader, const glm::mat4& transform) const {
    for (size_t i = 0; i < objects.size(); i++) {
        if (objects.at(i).name == "bond") meshes.at("bond").render(shader, transform * objects.at(i).getModel());
        else meshes.at(objects.at(i).name).render(shader, transform * objects.at(i).getModel());
    }
}

/*
Returns the number of objects (atoms and bonds).
*/
size_t Geometry::size() const {
    return objects.size();
}
