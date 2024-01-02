#!/usr/bin/env python

import numpy as np

import sys

if __name__ == "__main__":
    # check the number of arguments
    if len(sys.argv) != 5: print("INVALID NUMBER OF ARGUMENTS"); sys.exit(1)

    # replace the exponential notation
    sys.argv[1] = sys.argv[1].replace("^", "**")

    # define the x values
    x = np.linspace(float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))

    # define the potential
    y = eval(sys.argv[1])

    # print the header
    print("#         x                     V")

    # print the potential
    for i in range(len(x)): print("{0:21.14f} {1:21.14f}".format(x[i], y[i]))
