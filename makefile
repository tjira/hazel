INCLUDE := -isystem lib/boost/install/include -isystem lib/eigen -isystem lib/glad/include -isystem lib/glfw/install/include -isystem lib/glm -isystem lib/imgui -isystem lib/libint/install/include
FLAGS := -std=c++17

ifeq ($(DEBUG), 1)
FLAGS += -flarge-source-files -g -MMD -MP -O0 -pedantic -Wall -Wextra
else
FLAGS += -fopenmp -O2
endif
FLAGS += -DIMGUI_DEFINE_MATH_OPERATORS -DGPPFLAGS="$(FLAGS)"

all: .build bin bin/hazel bin/hview
libs: boost eigen glad glfw glm imgui libint

# Include Dependencies =================================================================================================

-include $(wildcard .build/*.d)

# Link =================================================================================================================

bin/hazel: .build/hazel.o .build/hartreefock.o .build/molecule.o .build/ptable.o .build/timer.o
	g++ $(FLAGS) $(INCLUDE) -o $@ $^ lib/boost/install/lib/libboost_program_options.a lib/libint/install/lib/libint2.a

bin/hview: .build/hview.o .build/buffer.o .build/gui.o .build/mesh.o .build/shader.o .build/glad.o .build/imgui.o .build/imgui_demo.o .build/imgui_dilog.o .build/imgui_draw.o .build/imgui_glfw.o .build/imgui_opengl.o .build/imgui_tables.o .build/imgui_widgets.o
	g++ $(FLAGS) $(INCLUDE) -o $@ $^ lib/boost/install/lib/libboost_program_options.a lib/glfw/install/lib/libglfw3.a -ldl

.build/hazel.o: hazel.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/hview.o: hview.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

# Project ==============================================================================================================

.build/buffer.o: src/buffer.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/gui.o: src/gui.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/hartreefock.o: src/hartreefock.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/mesh.o: src/mesh.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/molecule.o: src/molecule.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/ptable.o: src/ptable.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/shader.o: src/shader.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/timer.o: src/timer.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

# Libraries ============================================================================================================

boost:
	git clone --recursive --depth 1 https://github.com/boostorg/boost.git lib/boost
	cd lib/boost && ./bootstrap.sh --prefix="$$PWD/install" --with-libraries=json,program_options && cd -
	cd lib/boost && ./b2 install && cd -

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

libint:
	git clone --depth 1 https://github.com/evaleev/libint.git lib/libint
	cd lib/libint && ./autogen.sh && cd -
	cd lib/libint && ./configure CPPFLAGS="-I$$PWD/../eigen" --prefix="$$PWD/install" --with-boost="$$PWD/../boost/install" --with-cxxgen-optflags="-O3" && cd -
	cd lib/libint && make && make install && cd -

.build/glad.o: lib/glad/src/gl.c
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/imgui.o: lib/imgui/imgui.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/imgui_demo.o: lib/imgui/imgui_demo.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/imgui_dilog.o: lib/imgui/dialog/ImGuiFileDialog.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/imgui_draw.o: lib/imgui/imgui_draw.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/imgui_glfw.o: lib/imgui/imgui_impl_glfw.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/imgui_opengl.o: lib/imgui/imgui_impl_opengl3.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/imgui_tables.o: lib/imgui/imgui_tables.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/imgui_widgets.o: lib/imgui/imgui_widgets.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

# Miscellaneous ========================================================================================================

.build bin:
	mkdir -p $@

clean:
	rm -rf .build .cache .clangd .makefile .vscode bin lib
