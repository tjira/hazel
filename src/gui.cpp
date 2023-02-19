#include "../include/gui.h"

Gui::Gui(GLFWwindow* window) : window(window) {
    ImGui::CreateContext();
    ImGui_ImplOpenGL3_Init("#version 420");
    ImGui_ImplGlfw_InitForOpenGL(this->window, true);
    ImGui::GetIO().IniFilename = nullptr;
}

Gui::~Gui() {
    ImGui_ImplGlfw_Shutdown();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui::DestroyContext();
}

void Gui::render(Movie& scene) {
    GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window);

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    if (pointer->flags.renderOptions) {
        ImGui::Begin("Render Options", &pointer->flags.renderOptions, ImGuiWindowFlags_AlwaysAutoResize);
        ImGui::SliderFloat("Ambient", &pointer->light.ambient, 0, 1);
        ImGui::SliderFloat("Diffuse", &pointer->light.diffuse, 0, 1);
        ImGui::SliderFloat("Specular", &pointer->light.specular, 0, 1);
        ImGui::SliderFloat("Shininess", &pointer->light.shininess, 1, 128);
        ImGui::End();
    }

    if (pointer->flags.info) {
        ImGui::SetNextWindowPos({ 0, 0 }); ImGui::Begin("info", &pointer->flags.info,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_AlwaysAutoResize |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoBringToFrontOnFocus |
            ImGuiWindowFlags_NoFocusOnAppearing
        );
        ImGui::Text("%.1f", ImGui::GetIO().Framerate);
        ImGui::End();
    }

    if (ImGuiFileDialog::Instance()->Display("Import Molecule", ImGuiWindowFlags_NoCollapse, { 512, 288 })) {
        if (ImGuiFileDialog::Instance()->IsOk()) {
            scene = Movie::LoadTrajectory(ImGuiFileDialog::Instance()->GetFilePathName());
        }
        ImGuiFileDialog::Instance()->Close();
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
