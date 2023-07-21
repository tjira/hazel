#!/bin/bash

# extract the matrix lines and it's size
MATRIX=$(sed -n '/OVERLAP/,/KINETIC/p' | head -n -2 | tail -n +2)
SIZE=$(echo "$MATRIX" | head -n 1 | awk -F'[^0-9]+' '{print $2}')

# remove previous splits and split the matrix
rm -f .S* && echo "$MATRIX" | split -l $(( $SIZE + 1 )) - .S

# define the AWK command and paste the splits
AWK='{for(C=1;C<=NF;C++)printf"% 19.14f ",$C;print""}'
paste .S* | tail -n +2 | awk "$AWK" && rm -f .S*
