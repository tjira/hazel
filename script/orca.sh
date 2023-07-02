#!/bin/bash

cleanup() {
    rm -f -- *.0 *.tmp orca.densities orca.engrad orca.gbw orca.hess orca.inp orca_property.txt
}

inputOrca() {
    [[ "$4" == "cisd" ]] && echo -e "! $4 $5 HCORE NOFROZENCORE SCFCONV8 ENGRAD NUMGRAD\n*xyzfile $2 $3 $1" && return
    echo -e "! $4 $5 HCORE NOFROZENCORE SCFCONV8 ENGRAD\n*xyzfile $2 $3 $1"
}

# extract the number of atoms
ATOMS=$(echo "$(wc -l < "$1") - 2" | bc)

# cd to the input directory and run orca
cd "$(dirname "$1")" && inputOrca "$(basename "$1")" "$2" "$3" "$4" "$5" > orca.inp && orca orca.inp > orca.out && cleanup

# extract gradient and energy
[[ "$4" == "cisd" || "$4" == "hf" ]] && GRADIENT=$(grep -A $((ATOMS + 2)) "GRADIENT" orca.out | tail -n "$ATOMS" | awk '{printf("% 16.12f % 16.12f % 16.12f\n", $4, $5, $6)}')
[[ "$4" == "mp2" ]] && GRADIENT=$(grep -A "$ATOMS" "MP2 gradient" orca.out | tail -n "$ATOMS" | awk '{printf("% 16.12f % 16.12f % 16.12f\n", $2, $3, $4)}')
ENERGY=$(grep "FINAL SINGLE POINT ENERGY" orca.out | awk '{print $5}')

# print results
echo "ENERGY: $ENERGY" && echo -e "GRADIENT:\n$GRADIENT"