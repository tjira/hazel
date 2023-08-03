#!/bin/bash

DATA=$(cat)

cat > plot.gp <<- EOM
stats "-" nooutput
$DATA
e

plot "-" with lines lc rgb "#4E79A7" lw 2 notitle; pause mouse close
$DATA
EOM

gnuplot plot.gp
