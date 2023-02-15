#!/bin/bash

BASIS="STO-3G"
METHOD="RHF"

rm -rf out && cp xyz/* .

for MOL in *.xyz; do
    echo -e "! $METHOD $BASIS HCORE LARGEPRINT\n*xyzfile 0 1 $MOL" > "${MOL%.*}.inp" && orca "${MOL%.*}.inp" | tee "${MOL%.*}.out"
done

mkdir -p out && mv -- *.out out && rm -rf -- *.densities *.engrad *.gbw *.inp *.opt *.tmp *.txt *.xyz
