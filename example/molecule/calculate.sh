#!/bin/bash

BASIS="STO-3G"
METHOD="RHF"

PRINT="
Print[P_OneElec] 1
Print[P_Overlap] 1
Print[P_KinEn] 1
Print[P_Density] 1
"

cp xyz/* .
for MOL in *.xyz; do
    echo -e "! $METHOD $BASIS HCORE ENGRAD\n%output${PRINT}end\n*xyzfile 0 1 $MOL" > "${MOL%.*}.inp" && orca "${MOL%.*}.inp" | tee "${MOL%.*}.out"
done
mkdir -p out && mv -- *.out out && rm -rf -- *.0 *.densities *.Energies *.engrad *.gbw *.hess *.hostnames *.inp *.lastint *.lastscf *.opt *.tmp *.txt *.xyz
