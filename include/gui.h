#include <dialog/ImGuiFileDialog.h>
#include <imgui_impl_opengl3.h>
#include <imgui_impl_glfw.h>
#include "mesh.h"

class Gui {
public:
    Gui(GLFWwindow* window); ~Gui();
    void render(Mesh& mesh);

private:
    GLFWwindow* window;
};
