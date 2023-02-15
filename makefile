INCLUDE := -isystem lib/boost/install/include -isystem lib/glm -isystem lib/libint/install/include
FLAGS := -std=c++17

ifeq ($(DEBUG), 1)
FLAGS += -flarge-source-files -g -MMD -MP -O0 -pedantic -Wall -Wextra
else
FLAGS += -fopenmp -O2
endif
FLAGS += -DGPPFLAGS="$(FLAGS)"

all: .build bin bin/hazel bin/hview
libs: boost eigen glm libint

# Include Dependencies =================================================================================================

-include $(wildcard .build/*.d)

# Link =================================================================================================================

bin/hazel: .build/hazel.o .build/hartreefock.o .build/molecule.o .build/ptable.o .build/timer.o
	g++ $(FLAGS) $(INCLUDE) -o $@ $^ lib/boost/install/lib/libboost_program_options.a lib/libint/install/lib/libint2.a

bin/hview: .build/hview.o .build/buffer.o .build/gui.o .build/mesh3D.o .build/shader.o .build/glad.o .build/imgui.o .build/imgui_demo.o .build/imgui_dilog.o .build/imgui_draw.o .build/imgui_glfw.o .build/imgui_opengl.o .build/imgui_tables.o .build/imgui_widgets.o
	g++ $(FLAGS) $(INCLUDE) -o $@ $^ lib/boost/install/lib/libboost_program_options.a -ldl -lglfw

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

.build/mesh3D.o: src/mesh3D.cpp
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

glm:
	git clone --depth 1 https://github.com/g-truc/glm lib/glm

libint:
	git clone --depth 1 https://github.com/evaleev/libint.git lib/libint
	cd lib/libint && ./autogen.sh && cd -
	cd lib/libint && ./configure CPPFLAGS="-I$$PWD/../eigen" --prefix="$$PWD/install" --with-boost="$$PWD/../boost/install" --with-cxxgen-optflags="-O3" && cd -
	cd lib/libint && make && make install && cd -

.build/glad.o: lib/glad/glad.c
	g++ -c -o $@ $<

.build/imgui.o: lib/imgui/imgui.cpp
	g++ -c -o $@ $<

.build/imgui_demo.o: lib/imgui/imgui_demo.cpp
	g++ -c -o $@ $<

.build/imgui_dilog.o: lib/imgui/ImGuiFileDialog.cpp
	g++ -c -o $@ $<

.build/imgui_draw.o: lib/imgui/imgui_draw.cpp
	g++ -c -o $@ $<

.build/imgui_glfw.o: lib/imgui/imgui_impl_glfw.cpp
	g++ -c -o $@ $<

.build/imgui_opengl.o: lib/imgui/imgui_impl_opengl3.cpp
	g++ -c -o $@ $<

.build/imgui_tables.o: lib/imgui/imgui_tables.cpp
	g++ -c -o $@ $<

.build/imgui_widgets.o: lib/imgui/imgui_widgets.cpp
	g++ -c -o $@ $<

# Miscellaneous ========================================================================================================

.build bin:
	mkdir -p $@

clean:
	rm -rf hazel .build .cache .clangd .makefile .vscode bin *.exe
