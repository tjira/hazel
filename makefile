FLAGS := -std=c++17

INCLUDE := -isystem lib/libint/install/include
LIBS := lib/libint/install/lib/libint2.a

ifeq ($(DEBUG),1)
FLAGS += -g -MMD -O0 -Wall -Wextra
else
FLAGS += -O2
endif

# -include $(wildcard .build/*.d)

all: .build hazel

# Link =================================================================================================================

hazel: .build/hazel.o .build/hartreefock.o .build/molecule.o .build/printer.o
	g++ $(FLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

.build/hazel.o: hazel.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

# Project ==============================================================================================================

.build/hartreefock.o: src/hartreefock.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/molecule.o: src/molecule.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/printer.o: src/printer.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

# Libraries

libint:
	git clone --depth 1 https://github.com/evaleev/libint.git lib/libint
	cd lib/libint && ./autogen.sh && ./configure --prefix="$$PWD/install" --with-cxxgen-optflags="-O3" && make && make install && cd -

# Miscellaneous ========================================================================================================

clean:
	rm -rf hazel .build .cache .clangd .makefile .vscode *.exe *.json

.build:
	mkdir -p .build
