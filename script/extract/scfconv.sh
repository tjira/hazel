#!/bin/bash

# extract the data
DATA=$(sed -n '/FOCK METHOD/,/FINAL/p' | head -n -2 | tail -n +8)
echo "$DATA" | awk '{printf("%2d % 20.14f\n", $1, $2);}'
