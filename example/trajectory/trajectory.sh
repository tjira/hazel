#!/bin/bash

BASIS="STO-3G"
METHOD="RHF"

TIMESTEP="0.5_fs"
TEMP="298_k"
STEPS=1000

mkdir -p temp && cp ../molecule/*.xyz temp

cd temp && for MOL in *.xyz; do
    echo -e "! $METHOD $BASIS HCORE MD\n%md timestep $TIMESTEP\ninitvel $TEMP\nthermostat NHC $TEMP\ndump position stride 1 filename \"${MOL%.*}_trj.xyz\"\nrun $STEPS\nend\n*xyzfile 0 1 $MOL" > "${MOL%.*}.inp" && orca "${MOL%.*}.inp" | tee "${MOL%.*}.out"
    while read -r LINE; do {
        echo "$LINE" | awk 'NF==4 {printf("%2s % 2.10f % 2.10f % 2.10f\n", $1, $2, $3, $4)}';
        echo "$LINE" | awk 'NF==13 {printf("%s %i\n", tolower($5), $6)}';
        echo "$LINE" | awk 'NF==1 {printf("%s\n", $1)}';} >> "${MOL%.*}_trj.xyz.tmp"
    done < "${MOL%.*}_trj.xyz" && mv "${MOL%.*}_trj.xyz.tmp" "../$MOL"
done && cd .. && rm -r temp
