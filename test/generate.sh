#!/bin/bash

SYSTEMS=("acetone" "ammonia" "ethane" "ethanol" "ethylene" "methane" "pyrrole" "water")
BASES=("mini" "3-21g" "sto-3g" "6-31g" "6-311g" "cc-pvdz")

for SYSTEM in ${SYSTEMS[@]}; do
    for BASIS in ${BASES[@]}; do
        # echo the source file name
        FILE="energy_${SYSTEM,,}_hf_${BASIS//-/}.cpp"; echo "$FILE"

        # create the source file
        cat template/energy.in > "$FILE"

        # calculate the expected energy
        ENERGY=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" hf -d 3 5 -m 1000 -t 1e-8 | grep FINAL | awk '{print $4}')

        # replace values in the source file
        sed -i "s/BASIS/${BASIS//-/}/g" $FILE && sed -i "s/\"${BASIS//-/}\"/\"${BASIS^^}\"/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/EVALUE/$ENERGY/g" "$FILE"
    done
done

rm -f gmon.out
