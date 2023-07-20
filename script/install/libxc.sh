#!/bin/bash

# clone the repository
git clone --depth 1 https://github.com/ElectronicStructureLibrary/libxc.git libxc

# compile libint
cd libxc && autoreconf -i && ./configure --prefix="$PWD/install" && make && make install && cd ..
