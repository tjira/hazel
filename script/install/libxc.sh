#!/bin/bash

# clone the repository
git clone --depth 1 https://gitlab.com/libxc/libxc.git libxc

# compile libint
cd libxc && autoreconf -i && ./configure --prefix="$PWD/install" && make -j2 && make install && cd ..
