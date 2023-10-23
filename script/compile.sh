#!/bin/bash

g++ -s -O2 -mavx -DCXXFLAGS="\"-s -O2 -mavx\"" -DDATADIR="\"$PWD\"" -Iinclude -Ilib -isystem/usr/include/eigen3 src/* -lint2 -o hazel
