#!/bin/bash

# clone the repository
git clone https://github.com/evaleev/libint.git libint

# compile libint
cd libint && ./autogen.sh && ./configure --prefix="$PWD/install" --enable-1body=1 --enable-eri=1 --with-cxxgen-optflags="-O2" && make && make install && cd ..