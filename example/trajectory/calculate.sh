#!/bin/bash

BASIS="STO-3G"
METHOD="RHF"

TIMESTEP="0.5_fs"
TEMP="298_k"
STEPS=100

cp ../molecule/xyz/* .
for MOL in *.xyz; do
    echo -e "! $METHOD $BASIS HCORE LARGEPRINT MD\n%md timestep $TIMESTEP\ninitvel $TEMP\nthermostat NHC $TEMP\ndump position stride 1 filename \"${MOL%.*}_trj.xyz\"\nrun $STEPS\nend\n*xyzfile 0 1 $MOL" > "${MOL%.*}.inp" && orca "${MOL%.*}.inp" | tee "${MOL%.*}.out"
done
mkdir -p out xyz && mv -- *.out out && mv -- *_trj.xyz xyz && rm -rf -- *.csv *.densities *.engrad *.gbw *.inp *.log *.mdrestart *.opt *.tmp *.txt *.xyz
