#include <dialog/ImGuiFileDialog.h>
#include <GLFW/glfw3.h>
#include <imgui_impl_opengl3.h>
#include <imgui_impl_glfw.h>
#include "glfwpointer.h"
#include "movie.h"

class Gui {
public:

    // Constructors and destructors
    Gui(GLFWwindow* window); ~Gui();

    // State functions
    void render(Movie& movie);

private:
    GLFWwindow* window;
};
