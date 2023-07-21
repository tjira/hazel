#!/bin/bash

sed '/#include ".*"/d ; /#pragma once/d' \
    lib/argparse.hpp                     \
    include/eigen.h                      \
    include/timer.h                      \
    include/system.h                     \
    include/integral.h                   \
    include/transform.h                  \
    include/hf.h                         \
    include/mp.h                         \
    include/ci.h                         \
    include/gradient.h                   \
    include/hessian.h                    \
    include/optimizer.h                  \
    include/distributor.h                \
    src/eigen.cpp                        \
    src/timer.cpp                        \
    src/integral.cpp                     \
    src/transform.cpp                    \
    src/hf.cpp                           \
    src/mp.cpp                           \
    src/ci.cpp                           \
    src/gradient.cpp                     \
    src/hessian.cpp                      \
    src/optimizer.cpp                    \
    src/system.cpp                       \
    src/distributor.cpp                  \
    src/main.cpp                         \
| sed "s/^ *// ; /^$/d ; /^\/\//d" > hazel.cpp
