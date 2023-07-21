#!/bin/bash

# extract the matrix lines and it's size
MATRIX=$(sed -n '/NUCLEAR/,/COULOMB/p' | head -n -2 | tail -n +2)
SIZE=$(echo "$MATRIX" | head -n 1 | awk -F'[^0-9]+' '{print $2}')

# remove previous splits and split the matrix
rm -f .V* && echo "$MATRIX" | split -l $(( $SIZE + 1 )) - .V

# define the AWK command and paste the splits
AWK='{for(C=1;C<=NF;C++)printf"% 19.14f ",$C;print""}'
paste .V* | tail -n +2 | awk "$AWK" && rm -f .V*
