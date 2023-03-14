#include <backends/imgui_impl_opengl3.h>
#include <backends/imgui_impl_glfw.h>
#include <ImGuiFileDialog.h>
#include "trajectory.h"

class Gui {
public:

    // Constructors and destructors
    Gui(GLFWwindow* window); ~Gui();

    // State functions
    void render(Trajectory& movie);

private:
    GLFWwindow* window;
};
