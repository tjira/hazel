#!/bin/bash

SYSTEMS=("ammonia" "ethane" "ethylene" "formaldehyde" "methane" "water")
BASES=("mini" "3-21g" "sto-3g" "6-31g" "cc-pvdz")

CORES=64

# Hartree-Fock
for SYSTEM in ${SYSTEMS[@]}; do
    for BASIS in ${BASES[@]}; do
        # echo the source file name
        FILE="energy_${SYSTEM,,}_hf_${BASIS//-/}.cpp"; echo "$FILE"

        # create the source file
        cat template/energy_hf.in > "$FILE"

        # calculate the expected energy
        ENERGY=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES hf -d 3 5 -m 1000 -t 1e-8 | grep FINAL | awk '{print $4}')

        # replace values in the source file
        sed -i "s/BASIS/${BASIS//-/}/g" $FILE && sed -i "s/\"${BASIS//-/}\"/\"${BASIS^^}\"/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/EVALUE/$ENERGY/g" "$FILE"
    done
done

# Hartree-Fock analytical gradient
for SYSTEM in ${SYSTEMS[@]}; do
    for BASIS in ${BASES[@]}; do
        # echo the source file name
        FILE="grad_${SYSTEM,,}_hf_${BASIS//-/}.cpp"; echo "$FILE"

        # create the source file
        cat template/grad_hf.in > "$FILE"

        # calculate the expected energy
        GRAD=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES hf -d 3 5 -m 1000 -t 1e-8 -g | sed -n "/ANAL/,/GRAD/p" | tail -n +5 | head -n -2)
        GRAD=$(echo $GRAD | tr "\n" " " | tr -s " " | sed "s/\s*$//g" | sed "s/ /, /g")

        # replace values in the source file
        sed -i "s/BASIS/${BASIS//-/}/g" $FILE && sed -i "s/\"${BASIS//-/}\"/\"${BASIS^^}\"/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/GRADVAL/$GRAD/g" "$FILE"
    done
done

# MP2
for SYSTEM in ${SYSTEMS[@]}; do
    for BASIS in ${BASES[@]}; do
        # echo the source file name
        FILE="energy_${SYSTEM,,}_mp2_${BASIS//-/}.cpp"; echo "$FILE"

        # create the source file
        cat template/energy_mp2.in > "$FILE"

        # calculate the expected energy
        ENERGY=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES hf -d 3 5 -m 1000 -t 1e-8 mp2 | grep "FINAL MP2" | awk '{print $4}')

        # replace values in the source file
        sed -i "s/BASIS/${BASIS//-/}/g" $FILE && sed -i "s/\"${BASIS//-/}\"/\"${BASIS^^}\"/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/EVALUE/$ENERGY/g" "$FILE"
    done
done

rm -f gmon.out
