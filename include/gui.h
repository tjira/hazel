#include <dialog/ImGuiFileDialog.h>
#include <imgui_impl_opengl3.h>
#include <imgui_impl_glfw.h>
#include "scene.h"

class Gui {
public:
    Gui(GLFWwindow* window); ~Gui();
    void render(Scene& scene);

private:
    GLFWwindow* window;
};
