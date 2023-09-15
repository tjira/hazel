#!/usr/bin/env python

import matplotlib.animation as anm
import matplotlib.pyplot as plt

import argparse as ap
import numpy as np

import sys


input, data = sys.stdin.read(), []
if __name__ == "__main__":
    parser = ap.ArgumentParser(prog="Hazel Plotter", description="Plotting script for the Hazel package.")
    parser.add_argument("--mdenergy", action="store_true")
    parser.add_argument("--hfconv", action="store_true")
    parser.add_argument("--pes", action="store_true")
    parser.add_argument("--wfn", action="store_true")
    args = parser.parse_args()

    if args.hfconv:
        block = input[input.find("HARTREE-FOCK"):input.find("NUCLEAR REPULSION")]
        data = [float(line.split()[1]) for line in block.split("\n")[7:-2]]

    if args.mdenergy:
        block = input[input.find("MOLECULAR DYNAMICS"):input.find("EXECUTION TIME:")]
        data = [float(line.split()[1]) for line in block.split("\n")[5:-3]]

    if args.pes:
        block = input[input.find("ENERGY SCAN"):input.find("EXECUTION TIME:")]
        data = [float(line.split()[1]) for line in block.split("\n")[4:-3]]

    if args.wfn:
        for line in [line for line in input.split("\n") if line]:
            data.append([]) if line[0] == "#" else data[-1].append([float(entry) for entry in line.split()])

    if data:
        fig, ax = plt.subplots()

        if isinstance(data[0], float):
            data = [[float(i), data[i - 1]] for i in range(1, len(data) + 1)]
            
        if isinstance(data[0][0], float):
            for i in range(1, len(data[0])):
                plt.plot(np.array(data).T[0], np.array(data).T[i])

        elif isinstance(data[0][0][0], float):
            plots = [plt.plot(np.array(data[0]).T[0], np.array(data[0]).T[i])[0] for i in range(1, len(data[0][0]))]
            update = lambda i: [plots[j - 1].set_ydata(np.array(data[i]).T[j]) for j in range(1, len(data[0][0]))]
            ani = anm.FuncAnimation(fig, update, frames=len(data), interval=30)

        plt.tight_layout(); plt.show(); plt.close("all")
