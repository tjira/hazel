#!/bin/bash

./script/potential.py "0.5*x^2" -16 16 1024 > pes.dat && hazel qd -n 3 --imaginary && cat wavefunction.dat | ./script/plotwfn.py
