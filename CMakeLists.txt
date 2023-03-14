cmake_minimum_required(VERSION 3.25)
cmake_policy(SET CMP0135 NEW)
set(CMAKE_CXX_STANDARD 20)
include(FetchContent)
project(Hazel)

# set output flags
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -pedantic -Wall -Wextra -Wno-deprecated")
set(CMAKE_CXX_FLAGS_RELEASE "-s -O2 -Wno-deprecated")

# set glfw cmake flags
set(GLFW_BUILD_EXAMPLES OFF)
set(GLFW_BUILD_TESTS OFF)
set(GLFW_BUILD_DOCS OFF)

# download tar releases of some libraries
FetchContent_Declare(imdialog URL https://github.com/aiekick/ImGuiFileDialog/archive/refs/tags/v0.6.3.tar.gz)
FetchContent_Declare(libint URL https://github.com/evaleev/libint/releases/download/v2.7.2/libint-2.7.2-mpqc4.tgz)

# git clone some libraries
FetchContent_Declare(argparse GIT_REPOSITORY https://github.com/p-ranav/argparse.git GIT_TAG v2.9)
FetchContent_Declare(glad GIT_REPOSITORY https://github.com/Dav1dde/glad.git GIT_TAG v2.0.4 SOURCE_SUBDIR cmake)
FetchContent_Declare(glfw GIT_REPOSITORY https://github.com/glfw/glfw.git GIT_TAG 3.3.8)
FetchContent_Declare(imgui GIT_REPOSITORY https://github.com/ocornut/imgui.git GIT_TAG v1.88)
FetchContent_Declare(json GIT_REPOSITORY https://github.com/nlohmann/json.git GIT_TAG v3.11.2)

# fetch the dependencies
FetchContent_MakeAvailable(argparse glad glfw imdialog imgui json libint)

# generate glad library
glad_add_library(glad API gl:core=4.2)

# add imgui backends library
add_library(imbackend ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)
target_include_directories(imbackend PRIVATE ${imgui_SOURCE_DIR})

# add imgui library
add_library(imgui ${imgui_SOURCE_DIR}/imgui.cpp ${imgui_SOURCE_DIR}/imgui_draw.cpp ${imgui_SOURCE_DIR}/imgui_tables.cpp ${imgui_SOURCE_DIR}/imgui_widgets.cpp)
target_include_directories(imgui PRIVATE ${imgui_SOURCE_DIR})

# add imdialog library
add_library(imdialog ${imdialog_SOURCE_DIR}/ImGuiFileDialog.cpp)
target_include_directories(imdialog PRIVATE ${imgui_SOURCE_DIR})

# add hview executable
add_executable(hview hview.cpp src/buffer.cpp src/geometry.cpp src/gui.cpp src/mesh.cpp src/ptable.cpp src/trajectory.cpp src/shader.cpp)
target_include_directories(hview PRIVATE ${imdialog_SOURCE_DIR} ${imgui_SOURCE_DIR})
target_link_libraries(hview argparse glad glfw imbackend imdialog imgui)

# add hazel executable
add_executable(hazel hazel.cpp src/hartreefock.cpp src/logger.cpp src/ptable.cpp src/system.cpp src/timer.cpp)
target_include_directories(hazel PRIVATE ${json_SOURCE_DIR}/include)
target_link_libraries(hazel argparse libint2)
