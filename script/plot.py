#!/usr/bin/env python

import matplotlib.pyplot as plt

import argparse as ap
import numpy as np

import sys


input, data = sys.stdin.read(), []
if __name__ == "__main__":
    parser = ap.ArgumentParser(prog="Hazel Plotter", description="Plotting script for the Hazel package.")
    parser.add_argument("--mdenergy", action="store_true")
    parser.add_argument("--hfconv", action="store_true")
    args = parser.parse_args()

    print(input.strip())

    if args.hfconv:
        block = input[input.find("HARTREE-FOCK"):input.find("NUCLEAR REPULSION")]
        data = [float(line.split()[1]) for line in block.split("\n")[7:-2]]

    if args.mdenergy:
        block = input[input.find("MOLECULAR DYNAMICS"):input.find("EXECUTION TIME:")]
        data = [float(line.split()[1]) for line in block.split("\n")[5:-3]]

    if data:
        plt.plot(data)
        plt.tight_layout()
        plt.show()
