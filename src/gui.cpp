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

void Gui::render(TrajectoryGraphic& trajectory) {
    GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window);

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    if (pointer->flags.options) {
        static int subdivisions = SUBDIVISIONS, sectors = SECTORS;
        static bool smooth = SMOOTH;

        auto remeshCylinders = [](int sectors, bool smooth) {
            MoleculeGraphic::meshes.at("bond") = Mesh::Cylinder(sectors, smooth, "bond"); 
        };
        auto remeshSpheres = [](int subdivisions, bool smooth) {
            for (auto& [symbol, object] : ptable) {
                MoleculeGraphic::meshes.at(symbol) = Mesh::Icosphere(subdivisions, smooth, symbol);
                MoleculeGraphic::meshes.at(symbol).setColor(object.color);
            }
        };

        ImGui::Begin("Options", &pointer->flags.options, ImGuiWindowFlags_AlwaysAutoResize);
        if (ImGui::Checkbox("Smooth", &smooth)) {
            remeshSpheres(subdivisions, smooth);
            remeshCylinders(sectors, smooth);
        }
        if (ImGui::SliderInt("Sphere", &subdivisions, 0, 6)) remeshSpheres(subdivisions, smooth);
        if (ImGui::SliderInt("Cylinder", &sectors, 4, 128)) remeshCylinders(sectors, smooth);
        ImGui::SliderFloat("Ambient", &pointer->light.ambient, 0, 1);
        ImGui::SliderFloat("Diffuse", &pointer->light.diffuse, 0, 1);
        ImGui::SliderFloat("Specular", &pointer->light.specular, 0, 1);
        ImGui::SliderFloat("Shininess", &pointer->light.shininess, 1, 128);

        ImGui::SliderInt("Frame", &trajectory.getFrame(), 0, trajectory.size() - 1);
        if (ImGui::SliderFloat("Binding Factor", &pointer->bindingFactor, 0.001f, 0.05f)) {
            for (auto& molecule : trajectory.getGeoms()) molecule.rebind(pointer->bindingFactor);
        }
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

    if (pointer->flags.pause) {
        ImGui::SetNextWindowPos({ (float)pointer->width - 58, 0 });
        ImGui::Begin("pause", &pointer->flags.pause,
            ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_AlwaysAutoResize |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoBringToFrontOnFocus |
            ImGuiWindowFlags_NoFocusOnAppearing
        );
        ImGui::Text("PAUSED");
        ImGui::End();
    }

    if (ImGuiFileDialog::Instance()->Display("Import Molecule", ImGuiWindowFlags_NoCollapse, { 512, 288 })) {
        if (ImGuiFileDialog::Instance()->IsOk()) {
            trajectory = TrajectoryGraphic::Load(ImGuiFileDialog::Instance()->GetFilePathName());
        }
        ImGuiFileDialog::Instance()->Close();
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
