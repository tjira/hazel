#!/bin/bash

if [[ $# -eq 0 ]]; then
    echo "No molecule file provided."; exit 1
fi

if [[ $# -eq 1 ]]; then
    echo "No charge provided."; exit 1
fi

if [[ $# -eq 2 ]]; then
    echo "No multiplicity provided."; exit 1
fi

if [[ $# -eq 3 ]]; then
    echo "No method name provided."; exit 1
fi

if [[ $# -eq 4 ]]; then
    echo "No basis name provided."; exit 1
fi

cleanup() {
    rm -f -- *.0 *.tmp orca.densities orca.engrad orca.gbw orca.hess orca.inp orca_property.txt orca_*
}

inputOrca() {
    if [[ "$6" == "grad" ]]; then
        [[ "$4" == "ccsd" ]] && echo -e "! $4 $5 HCORE NOFROZENCORE SCFCONV8 ENGRAD NUMGRAD\n*xyzfile $2 $3 $1" && return
        [[ "$4" == "cisd" ]] && echo -e "! $4 $5 HCORE NOFROZENCORE SCFCONV8 ENGRAD NUMGRAD\n*xyzfile $2 $3 $1" && return
        echo -e "! $4 $5 HCORE NOFROZENCORE SCFCONV8 ENGRAD\n*xyzfile $2 $3 $1"
    else
        echo -e "! $4 $5 HCORE NOFROZENCORE SCFCONV8\n*xyzfile $2 $3 $1"
    fi
}

# extract the number of atoms
ATOMS=$(echo "$(wc -l < "$1") - 2" | bc)

# cd to the input directory and run orca
cd "$(dirname "$1")" && inputOrca "$(basename "$1")" "$2" "$3" "$4" "$5" "$6" > orca.inp && orca orca.inp > orca.out && cleanup

# extract gradient
[[ "$4" == "ccsd" ]] && GRADIENT=$(grep -A $((ATOMS + 1)) "GRADIENT" orca.out | tail -n "$ATOMS" | awk '{printf("% 16.12f % 16.12f % 16.12f\n", $4, $5, $6)}')
[[ "$4" == "cisd" ]] && GRADIENT=$(grep -A $((ATOMS + 1)) "GRADIENT" orca.out | tail -n "$ATOMS" | awk '{printf("% 16.12f % 16.12f % 16.12f\n", $4, $5, $6)}')
[[ "$4" == "hf" ]] && GRADIENT=$(grep -A $((ATOMS + 2)) "GRADIENT" orca.out | tail -n "$ATOMS" | awk '{printf("% 16.12f % 16.12f % 16.12f\n", $4, $5, $6)}')
[[ "$4" == "mp2" ]] && GRADIENT=$(grep -A "$ATOMS" "MP2 gradie" orca.out | tail -n "$ATOMS" | awk '{printf("% 16.12f % 16.12f % 16.12f\n", $2, $3, $4)}')
[[ "$4" == "fci" ]] && GRADIENT=$(grep -A $((ATOMS + 2)) "GRADIENT" orca.out | tail -n "$ATOMS" | awk '{printf("% 16.12f % 16.12f % 16.12f\n", $4, $5, $6)}')

# extract energy
[[ "$4" == "ccsd" ]] && ENERGY=$(grep "FINAL SINGLE POINT ENERGY" orca.out | awk '{print $5}')
[[ "$4" == "cisd" ]] && ENERGY=$(grep "FINAL SINGLE POINT ENERGY" orca.out | awk '{print $5}')
[[ "$4" == "hf" ]] && ENERGY=$(grep "FINAL SINGLE POINT ENERGY" orca.out | awk '{print $5}')
[[ "$4" == "mp2" ]] && ENERGY=$(grep "FINAL SINGLE POINT ENERGY" orca.out | awk '{print $5}')
[[ "$4" == "fci" ]] && ENERGY=$(grep "FINAL SINGLE POINT ENERGY" orca.out | tail -n 1 | awk '{print $5}')

# print results
echo "ENERGY: $ENERGY" && [[ "$6" == "grad" ]] && echo -e "GRADIENT:\n$GRADIENT"
