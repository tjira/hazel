#include "mesh3D.h"

#include <imgui_impl_opengl3.h>
#include <imgui_impl_glfw.h>
#include <dialog/ImGuiFileDialog.h>

class Gui {
public:
    Gui(GLFWwindow* window); ~Gui();
    void render(Mesh3D& mesh);

private:
    GLFWwindow* window;
};
