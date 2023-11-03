# set glfw cmake flags
set(GLFW_BUILD_EXAMPLES OFF)
set(GLFW_BUILD_TESTS OFF)
set(GLFW_BUILD_DOCS OFF)

# declare library sources libraries
FetchContent_Declare(imdialog SYSTEM GIT_REPOSITORY https://github.com/aiekick/ImGuiFileDialog.git GIT_TAG 8371c8612ddfe9b23464ba9b02770ea046e81f98)
FetchContent_Declare(implot SYSTEM GIT_REPOSITORY https://github.com/epezent/implot.git GIT_TAG 065acc3319f0422479c0fed5a5edccd0f563729f)
FetchContent_Declare(imgui SYSTEM GIT_REPOSITORY https://github.com/ocornut/imgui.git GIT_TAG 1ab63d925f21e03be7735661500e5b914dd93c19)
FetchContent_Declare(glad SYSTEM GIT_REPOSITORY https://github.com/Dav1dde/glad.git GIT_TAG 2348b07c1ab4504d60398713781d8a57880234fa)
FetchContent_Declare(stb SYSTEM GIT_REPOSITORY https://github.com/nothings/stb.git GIT_TAG beebb24b945efdea3b9bba23affb8eb3ba8982e7)
FetchContent_Declare(glfw SYSTEM GIT_REPOSITORY https://github.com/glfw/glfw.git GIT_TAG 3eaf1255b29fdf5c2895856c7be7d7185ef2b241)
FetchContent_Declare(glm SYSTEM GIT_REPOSITORY https://github.com/g-truc/glm.git GIT_TAG 47585fde0c49fa77a2bf2fb1d2ead06999fd4b6e)

# fetch the libraries
FetchContent_MakeAvailable(glad glfw glm imdialog imgui implot stb)

# generate glad library
add_subdirectory(${glad_SOURCE_DIR}/cmake)
glad_add_library(glad API gl:core=4.2)

# create an alternate directory for glad loading (mainly used by libraries)
add_custom_command(TARGET glad POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_BINARY_DIR}/gladsources/glad/include/glad/gl.h
    ${CMAKE_CURRENT_BINARY_DIR}/gladsources/glad/include/GL/gl.h
)

# include headers for imdialog library
target_include_directories(ImGuiFileDialog PUBLIC ${imgui_SOURCE_DIR})

# include libraries for the hview binary
include_directories(tool/hview/include ${argparse_SOURCE_DIR}/include ${imgui_SOURCE_DIR} ${implot_SOURCE_DIR} ${stb_SOURCE_DIR})

# add hview executable
add_executable(hview
    tool/hview/src/buffer.cpp
    tool/hview/src/geometry.cpp
    tool/hview/src/gui.cpp
    tool/hview/src/main.cpp
    tool/hview/src/mesh.cpp
    tool/hview/src/ptable.cpp
    tool/hview/src/shader.cpp
    tool/hview/src/trajectory.cpp

    # imgui backends
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp

    # imgui
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp

    # implot
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_demo.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp
)

# link hview executable
target_link_libraries(hview glad glfw glm::glm ImGuiFileDialog)
