#!/bin/bash

# extract the matrix lines and it's size
MATRIX=$(sed -n '/COULOMB/,/--/p' | head -n -2 | tail -n +2)
SIZE=$(echo "$MATRIX" | head -n 1 | awk -F'[^0-9]+' '{print $2}')

# remove previous splits and split the matrix
rm -f .J* && echo "$MATRIX" | split -l $(( $SIZE + 1 )) - .J

# define the AWK command and paste the splits
AWK='{for(C=1;C<=NF;C++)printf"% 20.14f ",$C;print""}'
paste .J* | tail -n +2 | awk "$AWK" && rm -f .J*
