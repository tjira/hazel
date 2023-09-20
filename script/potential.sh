#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo "No function provided."; exit 1
fi

if [[ $# -eq 1 ]]; then
    echo "No range provided."; exit 1
fi

if [[ $# -eq 2 ]]; then
    echo "No range provided."; exit 1
fi

if [[ $# -eq 3 ]]; then
    echo "No number of points provided."; exit 1
fi

echo "#     X           E"
for X in $(seq "$2" $(echo "($3 - $2) / ($4 - 1)" | bc -l) "$3"); do
    V=$(echo "${1/x/$X}" | bc -l)
    printf "%.8f %.8f\n" $X $V
done
