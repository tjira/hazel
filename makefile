FLAGS := -std=c++17

LIBS := lib/libint/install/lib/libint2.a -l boost_json -l boost_program_options
INCLUDE := -isystem lib/libint/install/include

ifeq ($(DEBUG), 1)
FLAGS += -flarge-source-files -g -MMD -MP -O0 -pedantic -Wall -Wextra
else
FLAGS += -fopenmp -O2
endif
FLAGS += -DGPPFLAGS="$(FLAGS)"

all: .build bin bin/hazel

# Include Dependencies =================================================================================================

-include $(wildcard .build/*.d)

# Link =================================================================================================================

bin/hazel: .build/hazel.o .build/hartreefock.o .build/molecule.o .build/ptable.o .build/timer.o
	g++ $(FLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

.build/hazel.o: hazel.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

# Project ==============================================================================================================

.build/hartreefock.o: src/hartreefock.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/molecule.o: src/molecule.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/ptable.o: src/ptable.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

.build/timer.o: src/timer.cpp
	g++ $(FLAGS) $(INCLUDE) -c -o $@ $<

# Libraries ============================================================================================================

libint:
	git clone --depth 1 https://github.com/evaleev/libint.git lib/libint && cd lib/libint && ./autogen.sh && cd -
	cd lib/libint && ./configure --prefix="$$PWD/install" --with-cxxgen-optflags="-O3" && cd -
	cd lib/libint && make && make install && cd -

# Miscellaneous ========================================================================================================

.build bin:
	mkdir -p $@

clean:
	rm -rf hazel .build .cache .clangd .makefile .vscode bin *.exe
