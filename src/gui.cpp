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

void Gui::render(Trajectory& trajectory) {
    GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window);

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    if (pointer->flags.options) {
        static float bindingFactor = BINDINGFACTOR, bondSize = BONDSIZE, atomSizeFactor = ATOMSIZEFACTOR;
        static int subdivisions = SUBDIVISIONS, sectors = SECTORS;
        static bool smooth = SMOOTH;

        auto remeshCylinders = [](int sectors, bool smooth) {
            Geometry::meshes.at("bond") = Mesh::Cylinder(sectors, smooth, "bond"); 
        };
        auto remeshSpheres = [](int subdivisions, bool smooth) {
            for (auto& [symbol, object] : ptable) {
                Geometry::meshes.at(symbol) = Mesh::Icosphere(subdivisions, smooth, symbol);
                Geometry::meshes.at(symbol).setColor(object.color);
            }
        };

        ImGui::Begin("Options", &pointer->flags.options, ImGuiWindowFlags_AlwaysAutoResize);

        if (ImGui::Checkbox("Smooth", &smooth)) {
            remeshSpheres(subdivisions, smooth);
            remeshCylinders(sectors, smooth);
        }

        ImGui::Separator();

        if (ImGui::SliderInt("Sphere", &subdivisions, 0, 6)) remeshSpheres(subdivisions, smooth);
        if (ImGui::SliderInt("Cylinder", &sectors, 4, 128)) remeshCylinders(sectors, smooth);
        if (ImGui::SliderFloat("Atom Size Factor", &atomSizeFactor, 0.001, 0.02)) {
            for (auto& molecule : trajectory.getGeoms()) molecule.setAtomSizeFactor(atomSizeFactor);
        }
        if (ImGui::SliderFloat("Bond Size", &bondSize, 0.01, 0.2)) {
            for (auto& molecule : trajectory.getGeoms()) molecule.setBondSize(bondSize);
        }

        ImGui::Separator();
        
        if (ImGui::SliderFloat("Binding Factor", &bindingFactor, 0, 0.05f)) {
            for (auto& molecule : trajectory.getGeoms()) molecule.rebind(bindingFactor);
        }

        ImGui::Separator();
        
        ImGui::SliderFloat("Ambient", &pointer->light.ambient, 0, 1);
        ImGui::SliderFloat("Diffuse", &pointer->light.diffuse, 0, 1);
        ImGui::SliderFloat("Specular", &pointer->light.specular, 0, 1);
        ImGui::SliderFloat("Shininess", &pointer->light.shininess, 1, 128);

        ImGui::Separator();

        ImGui::SliderInt("Frame", &trajectory.getFrame(), 0, trajectory.size() ? trajectory.size() - 1 : 0);
        ImGui::SliderFloat("Timeout", &trajectory.getWait(), 0.001, 16);

        ImGui::Separator();
        
        if (ImGui::Button("Center")) {
            trajectory.moveBy(-trajectory.getGeoms().at(trajectory.getFrame()).getCenter());
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
            trajectory = Trajectory::Load(ImGuiFileDialog::Instance()->GetFilePathName());
        }
        ImGuiFileDialog::Instance()->Close();
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
