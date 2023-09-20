#!/usr/bin/env python

import matplotlib.animation as anm
import matplotlib.pyplot as plt

import argparse as ap
import numpy as np

import sys
import re


input, data, pot, energy = sys.stdin.read(), [], [], []
xmin, xmax, ymin, ymax = 100, -100, 100, -100

if __name__ == "__main__":
    parser = ap.ArgumentParser(prog="Hazel Plotter", description="Plotting script for the Hazel package.")
    parser.add_argument("--mdenergy", action="store_true")
    parser.add_argument("--hfconv", action="store_true")
    parser.add_argument("--pes", action="store_true")
    parser.add_argument("--wfn", action="store_true")
    parser.add_argument("--pot", default="pes.dat")
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
        with open(args.pot, "r") as file:
            pot = np.array([[float(cell) for cell in line.split()] for line in file.read().strip().split("\n")[1:]])
        energy = np.array([float(re.sub("[E=,\[\]]", "", cell)) for cell in input.split("\n")[0].split()[len(data[0][0]) + 1:]])

    if data:
        fig, ax = plt.subplots()

        if isinstance(data[0], float):
            data = [[float(i), data[i - 1]] for i in range(1, len(data) + 1)]
            
        if args.hfconv or args.mdenergy or args.pes:
            for i in range(1, len(data[0])):
                plt.plot(np.array(data).T[0], np.array(data).T[i])

        if args.wfn:
            if len(pot): plt.plot(np.array(pot).T[0], np.array(pot).T[1])
            for i in range(len(data)):
                xming = np.min(np.array(data[0]).T[0][np.apply_along_axis(lambda row: row.sum(), 1, np.abs(np.array(data[i])[:, 1:])) > 0.01])
                xmaxg = np.max(np.array(data[0]).T[0][np.apply_along_axis(lambda row: row.sum(), 1, np.abs(np.array(data[i])[:, 1:])) > 0.01])
                yming = np.min(np.array(data[i])[:, 1:3]) + np.min(energy)
                ymaxg = np.max(np.array(data[i])[:, -2:]) + np.max(energy)
                xmin = xming if xming < xmin else xmin
                xmax = xmaxg if xmaxg > xmax else xmax
                ymin = yming if yming < ymin else ymin
                ymax = ymaxg if ymaxg > ymax else ymax
            ymin = np.min([ymin, np.min(pot[:, 1])])
            plots = [plt.plot(np.array(data[0]).T[0], np.array(data[0]).T[i] + energy[(i - 1) // 2])[0] for i in range(1, len(data[0][0]))]
            update = lambda i: [plots[j - 1].set_ydata(np.array(data[i]).T[j] + energy[(j - 1) // 2]) for j in range(1, len(data[0][0]))]
            plt.axis([xmin, xmax, ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)])
            ani = anm.FuncAnimation(fig, update, frames=len(data), interval=30);

        plt.tight_layout(); plt.show(); plt.close("all")
