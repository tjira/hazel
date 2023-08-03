#!/bin/bash

# extract the data
DATA=$(sed -n '/MOLECULAR DYN/,/FINAL/p' | head -n -4 | tail -n +6)
echo "$DATA" | awk '{printf("%2d % 20.14f\n", $1, $2);}'
