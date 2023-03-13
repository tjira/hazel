INCLUDE := lib/argparse/include lib/boost lib/eigen lib/glad/include lib/glfw/install/include lib/glm lib/imgui lib/json/single_include lib/libint/install/include
FLAGS := -std=c++20 -MMD -MP

# Variable Modifications ===============================================================================================

ifeq ($(DEBUG), 1)
FLAGS += -flarge-source-files -g -O0 -pedantic -Wall -Wextra
else
FLAGS += -fopenmp -O2
endif
FLAGS += -DIMGUI_DEFINE_MATH_OPERATORS -DGPPFLAGS="$(FLAGS)"
INCLUDE := $(addprefix -isystem , $(INCLUDE))

# Functions =============================================================================================================

uniq = $(if $1,$(firstword $1) $(call uniq,$(filter-out $(firstword $1),$1)))

# Object Files ==========================================================================================================

IMGUI := imgui.o imgui_demo.o imgui_dilog.o imgui_draw.o imgui_glfw.o imgui_opengl.o imgui_tables.o imgui_widgets.o
HAZEL := hartreefock.o logger.o ptable.o system.o timer.o
HVIEW := buffer.o geometry.o gui.o mesh.o ptable.o shader.o trajectory.o

# Targets ===============================================================================================================

all: bin build bin/hazel bin/hview
libs: argparse boost eigen glad glfw glm imgui json libint

# Include Dependencies ==================================================================================================

-include $(wildcard build/*.d)

# Link ==================================================================================================================

bin/hazel: $(addprefix build/, hazel.o $(HAZEL))
	g++ $(FLAGS) $(INCLUDE) -o $@ $^ lib/libint/install/lib/libint2.a

bin/hview: $(addprefix build/, hview.o glad.o $(HVIEW) $(IMGUI))
	g++ $(FLAGS) $(INCLUDE) -o $@ $^ lib/glfw/install/lib/libglfw3.a -ldl

# Project Files =========================================================================================================

build/hazel.o build/hview.o: build/%.o: %.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

$(call uniq, $(addprefix build/, $(HAZEL) $(HVIEW))): build/%.o: src/%.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

# Libraries =============================================================================================================

argparse:
	git clone --depth 1 https://github.com/p-ranav/argparse.git lib/argparse

boost:
	git clone --recursive --depth 1 https://github.com/boostorg/boost.git lib/boost
	cd lib/boost && ./bootstrap.sh && ./b2 headers && cd -

eigen:
	git clone --depth 1 https://gitlab.com/libeigen/eigen.git lib/eigen

glad:
	wget -q "https://gen.glad.sh"$$(curl https://gen.glad.sh/generate -s -X POST --data-raw 'generator=c&api=gl%3D4.2&profile=gl%3Dcore&options=LOADER' | grep -o "\".*\"" | tail -n 1 | tr -d '""')"/glad.zip"
	mkdir -p lib/glad && unzip glad.zip -d lib/glad && rm glad.zip && cp -r lib/glad/include/glad lib/glad/include/GL

glfw:
	git clone --depth 1 https://github.com/glfw/glfw.git lib/glfw
	cd lib/glfw && cmake -B build -DCMAKE_INSTALL_PREFIX="$$PWD/install" -DGLFW_BUILD_DOCS=OFF -DGLFW_BUILD_EXAMPLES=OFF -DGLFW_BUILD_TESTS=OFF && cd -
	cd lib/glfw/build && make && make install && cd -

glm:
	git clone --depth 1 https://github.com/g-truc/glm.git lib/glm

imgui:
	git clone --depth 1 https://github.com/ocornut/imgui.git lib/imgui
	git clone --branch Lib_Only --depth 1 https://github.com/aiekick/ImGuiFileDialog.git lib/imgui/dialog
	cp lib/imgui/backends/imgui_impl_glfw* lib/imgui/backends/imgui_impl_opengl3* lib/imgui/

json:
	git clone --depth 1 https://github.com/nlohmann/json.git lib/json

libint:
	git clone --depth 1 https://github.com/evaleev/libint.git lib/libint
	cd lib/libint && ./autogen.sh && sed -i '63d;65,67d' include/libint2/cgshell_ordering.h && cd -
	cd lib/libint && ./configure CXX=g++ CPPFLAGS="-I$$PWD/../eigen" --prefix="$$PWD/install" --with-boost="$$PWD/../boost/install" --with-cxxgen-optflags="-O3" --with-cartgauss-ordering=orca && cd -
	cd lib/libint && make && make install && cd -

build/glad.o: lib/glad/src/gl.c
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

build/imgui.o: lib/imgui/imgui.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

build/imgui_demo.o: lib/imgui/imgui_demo.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

build/imgui_dilog.o: lib/imgui/dialog/ImGuiFileDialog.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

build/imgui_draw.o: lib/imgui/imgui_draw.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

build/imgui_glfw.o: lib/imgui/imgui_impl_glfw.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

build/imgui_opengl.o: lib/imgui/imgui_impl_opengl3.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

build/imgui_tables.o: lib/imgui/imgui_tables.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

build/imgui_widgets.o: lib/imgui/imgui_widgets.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

# Miscellaneous =========================================================================================================

bin build:
	@mkdir -p $@

clean:
	rm -rf build .cache .clangd .makefile .vscode bin lib
