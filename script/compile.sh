#!/bin/bash

g++ -s -O2 -DCXXFLAGS="\"-s O2\"" -DDATADIR="\"$PWD\"" -Iinclude -Ilib -isystem/usr/include/eigen3 src/* -lint2 -lxc -o hazel
