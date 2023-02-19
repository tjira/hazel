#include "include/gui.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

std::string vertex = R"(
#version 420 core
layout(location = 0) in vec3 i_position;
layout(location = 1) in vec3 i_normal;
layout(location = 2) in vec3 i_color;
uniform mat4 u_model, u_view, u_proj;
out vec3 fragment, normal, color;
out mat3 transform;
void main() {
    normal = normalize(mat3(transpose(inverse(u_model))) * i_normal);
    fragment = vec3(u_model * vec4(i_position, 1)), color = i_color;
    gl_Position = u_proj * u_view * vec4(fragment, 1);
    transform = inverse(mat3(u_view));
})";

std::string fragment = R"(
#version 420 core
struct Light { vec3 position; float ambient, diffuse, specular, shininess; };
uniform Light u_light; uniform vec3 u_camera;
in vec3 fragment, normal, color;
in mat3 transform;
out vec4 o_color;
void main() {
    vec3 lightPos = transform * u_light.position;
    vec3 reflection = reflect(-normalize(lightPos), normal);
    vec3 direction = normalize(u_camera - fragment);
    vec3 specular = vec3(pow(max(dot(direction, reflection), 0), u_light.shininess));
    vec3 diffuse = vec3(max(dot(normal, normalize(lightPos)), 0));
    o_color = vec4((vec3(u_light.ambient) + u_light.diffuse * diffuse + u_light.specular * specular), 1) * vec4(color, 1);
})";

void keyCallback(GLFWwindow* window, int key, int, int action, int mods) {
    if (GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window); action == GLFW_PRESS) {
        if (mods == GLFW_MOD_CONTROL) {
            if (key == GLFW_KEY_O) {
                std::string files = "Molecule Files{.xyz},All Files{.*}";
                ImGuiFileDialog::Instance()->OpenDialog("Import Molecule", "Import Molecule", files.c_str(), ".");
            } else if (key == GLFW_KEY_Q) {
                glfwSetWindowShouldClose(window, GLFW_TRUE);
            }
        }
        else if (key == GLFW_KEY_F1) pointer->flags.renderOptions = !pointer->flags.renderOptions;
        else if (key == GLFW_KEY_F11) {
            static int xpos0, ypos0, width0, height0;
            int xpos, ypos, width, height;
            if (pointer->flags.fullscreen = !pointer->flags.fullscreen; pointer->flags.fullscreen) {
                glfwGetWindowSize(pointer->window, &width0, &height0);
                glfwGetWindowPos(pointer->window, &xpos0, &ypos0);
                glfwGetMonitorWorkarea(glfwGetPrimaryMonitor(), &xpos, &ypos, &width, &height);
                glfwSetWindowMonitor(pointer->window, glfwGetPrimaryMonitor() , 0, 0, width, height, GLFW_DONT_CARE);
            } else {
                glfwSetWindowMonitor(pointer->window, nullptr , xpos0, ypos0, width0, height0, GLFW_DONT_CARE);
            }
        }
        else if (key == GLFW_KEY_F12) pointer->flags.info = !pointer->flags.info;
    }
}

void positionCallback(GLFWwindow* window, double x, double y) {
    GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window);
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS && !ImGui::GetIO().WantCaptureMouse) {
        glm::vec3 xaxis = glm::inverse(glm::mat3(pointer->camera.view)) * glm::vec3(0, 1, 0);
        glm::vec3 yaxis = glm::inverse(glm::mat3(pointer->camera.view)) * glm::vec3(1, 0, 0);
        pointer->camera.view = glm::rotate(pointer->camera.view, 0.01f * ((float)y - pointer->mouse.y), yaxis);
        pointer->camera.view = glm::rotate(pointer->camera.view, 0.01f * ((float)x - pointer->mouse.x), xaxis);
    }
    pointer->mouse = { x, y };
}

void resizeCallback(GLFWwindow* window, int width, int height) {
    if (GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window); width > 0 && height > 0) {
        pointer->camera.proj = glm::perspective(glm::radians(45.0f), (float)width / height, 0.01f, 1000.0f);
        pointer->width = width, pointer->height = height; glViewport(0, 0, width, height);
    }
}

void scrollCallback(GLFWwindow* window, double, double dy) {
    if (!ImGui::GetIO().WantCaptureMouse) {
        ((GLFWPointer*)glfwGetWindowUserPointer(window))->camera.view *= glm::mat4(glm::mat3(1.0f + 0.08f * (float)dy));
    }
}

void set(const Shader& shader, const GLFWPointer::Camera& camera, const GLFWPointer::Light& light) {
    shader.set<glm::vec3>("u_camera", -glm::inverse(glm::mat3(camera.view)) * glm::vec3(camera.view[3]));
    shader.set<glm::mat4>("u_view", camera.view);
    shader.set<glm::mat4>("u_proj", camera.proj);
    shader.set<glm::vec3>("u_light.position", light.position);
    shader.set<float>("u_light.shininess", light.shininess);
    shader.set<float>("u_light.specular", light.specular);
    shader.set<float>("u_light.ambient", light.ambient);
    shader.set<float>("u_light.diffuse", light.diffuse);
}

int main(int argc, char** argv) {
    // initialize the argument parser and container for the arguments
    po::options_description desc("options");
    po::positional_options_description pos;
    po::variables_map vm;

    // add options to the parser
    desc.add_options()
        ("input", po::value<std::string>(), "input file")
        ("version,v", "print version string")
        ("help,h", "produce help message")
    ;pos.add("input", 1);

    // extract the variables from the command line
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm); po::notify(vm);

    // print help if the help flag was provided
    if (vm.count("help")) {
        std::cout << desc << std::endl; return EXIT_SUCCESS;
    }

    // open the provided input
    if (!vm["input"].empty() && !std::filesystem::exists(vm["input"].as<std::string>())) {
        throw std::runtime_error("Input file does not exist.");
    }

    // Create GLFW variable struct
    GLFWPointer pointer; 

    // Initialize GLFW and throw error if failed
    if(!glfwInit()) {
        throw std::runtime_error("Error during GLFW initialization.");
    }

    // Pass OpenGL version and other hints
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, pointer.major);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, pointer.minor);
    glfwWindowHint(GLFW_SAMPLES, pointer.samples);

    // Create the window
    if (pointer.window = glfwCreateWindow(pointer.width, pointer.height, pointer.title.c_str(), nullptr, nullptr); !pointer.window) {
        throw std::runtime_error("Error during window creation.");
    }

    // Initialize GLAD
    if (glfwMakeContextCurrent(pointer.window); !gladLoadGL(glfwGetProcAddress)) {
        throw std::runtime_error("Error during GLAD initialization.");
    }

    // Enable some options
    glEnable(GL_DEPTH_TEST), glEnable(GL_CULL_FACE);
    glfwSetWindowUserPointer(pointer.window, &pointer);
    glfwSwapInterval(1);

    // Set event callbacks
    glfwSetCursorPosCallback(pointer.window, positionCallback);
    glfwSetWindowSizeCallback(pointer.window, resizeCallback);
    glfwSetScrollCallback(pointer.window, scrollCallback);
    glfwSetKeyCallback(pointer.window, keyCallback);

    // Initialize camera matrices
    pointer.camera.proj = glm::perspective(glm::radians(45.0f), (float)pointer.width / pointer.height, 0.01f, 1000.0f);
    pointer.camera.view = glm::lookAt({ 0.0f, 0.0f, 5.0f }, glm::vec3(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));

    {
        // Create scene, shader and GUI
        Movie movie;
        if (!vm["input"].empty()) {
            movie = Movie::LoadTrajectory(vm["input"].as<std::string>());
        }
        Shader shader(vertex, fragment);
        Gui gui(pointer.window);
        
        // Enter the render loop
        while (!glfwWindowShouldClose(pointer.window)) {
            
            // Clear the color and depth buffer
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // Set shade variables
            set(shader, pointer.camera, pointer.light);

            // Render mesh and GUI
            movie.render(shader);
            gui.render(movie);
            
            // Swap buffers and poll events
            glfwSwapBuffers(pointer.window);
            glfwPollEvents();
        }
    }

    // Clean up GLFW
    glfwTerminate();
}
