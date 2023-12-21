#!/usr/bin/env python

import matplotlib.animation as anm
import matplotlib.pyplot as plt

import argparse as ap
import numpy as np

import sys
import re

xmin, xmax, ymin, ymax = 100, -100, 100, -100
data, input = list(), sys.stdin.read()

if __name__ == "__main__":
    # create the argument parser
    parser = ap.ArgumentParser(prog="Hazel Wavefunction Plotter", description="Wavefunction plotting script for the Hazel package.")

    # add arguments
    parser.add_argument("--density", action="store_true")

    # parse arguments
    args = parser.parse_args()

    # parse the wavefunctions
    for line in [line for line in input.split("\n") if line]:
        data.append([]) if line[0] == "#" else data[-1].append([float(entry) for entry in line.split()])

    # load the potential file
    with open("pes.dat", "r") as file:
        pot = np.array([[float(cell) for cell in line.split()] for line in file.read().strip().split("\n")[1:]])

    # extract the energy from the wavefunction file
    energy = np.array([float(re.sub("[E=,\[\]]", "", cell)) for cell in input.split("\n")[0].split()[len(data[0][0]) + 1:]])

    # create the figure
    fig, ax = plt.subplots()

    # plot the potential
    plt.plot(np.array(pot).T[0], np.array(pot).T[1])

    # extract the wavefunction limits
    for i in range(len(data)):
        xming = np.min(np.array(data[0]).T[0][np.apply_along_axis(lambda row: row.sum(), 1, np.abs(np.array(data[i])[:, 1:])) > 0.01])
        xmaxg = np.max(np.array(data[0]).T[0][np.apply_along_axis(lambda row: row.sum(), 1, np.abs(np.array(data[i])[:, 1:])) > 0.01])
        yming = np.min(np.array(data[i])[:, 1:3]) + np.min(energy)
        ymaxg = np.max(np.array(data[i])[:, -2:]) + np.max(energy)
        xmin = xming if xming < xmin else xmin
        xmax = xmaxg if xmaxg > xmax else xmax
        ymin = yming if yming < ymin else ymin
        ymax = ymaxg if ymaxg > ymax else ymax

    # set the Y axis minimum
    ymin = np.min([ymin, np.min(pot[:, 1])])

    # set the plot limits
    plt.axis([xmin, xmax, ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)])

    # define the plots and animation update function
    plots = [plt.plot(np.array(data[0]).T[0], np.array(data[0]).T[i] + energy[(i - 1) // 2])[0] for i in range(1, len(data[0][0]))]
    update = lambda i: [plots[j - 1].set_ydata(np.array(data[i]).T[j] + energy[(j - 1) // 2]) for j in range(1, len(data[0][0]))]

    # create the animation
    ani = anm.FuncAnimation(fig, update, frames=len(data), interval=30);

    # show the plot
    plt.tight_layout(); plt.show(); plt.close("all")
