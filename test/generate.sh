#!/bin/bash

SYSTEMS=("ammonia" "ethane" "ethylene" "formaldehyde" "methane" "water")
BASES=("mini" "3-21g" "sto-3g" "6-31g" "6-31g*" "cc-pvdz")

CORES=64

# overlap matrix
for SYSTEM in "${SYSTEMS[@]}"; do
    for BASIS in "${BASES[@]}"; do
        # echo the source file name
        FILE="int_overlap_${SYSTEM,,}_$BASIS.cpp"; echo "$FILE"

        # create the source file
        cat template/int_overlap.in > "$FILE"

        # calculate the expected energy
        S=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES -p s | ../script/extract/overlap.sh)
        S=$(echo "$S" | tr "\n" " " | tr -s " " | sed "s/^ // ; s/\s*$//g ; s/ /, /g")

        # replace values in the source file
        sed -i "s/BNAME/$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/SVAL/$S/g" "$FILE"
    done
done

# kinetic matrix
for SYSTEM in "${SYSTEMS[@]}"; do
    for BASIS in "${BASES[@]}"; do
        # echo the source file name
        FILE="int_kinetic_${SYSTEM,,}_$BASIS.cpp"; echo "$FILE"

        # create the source file
        cat template/int_kinetic.in > "$FILE"

        # calculate the expected energy
        T=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES -p t | ../script/extract/kinetic.sh)
        T=$(echo "$T" | tr "\n" " " | tr -s " " | sed "s/^ // ; s/\s*$//g ; s/ /, /g")

        # replace values in the source file
        sed -i "s/BNAME/$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/TVAL/$T/g" "$FILE"
    done
done

# nuclear matrix
for SYSTEM in "${SYSTEMS[@]}"; do
    for BASIS in "${BASES[@]}"; do
        # echo the source file name
        FILE="int_nuclear_${SYSTEM,,}_$BASIS.cpp"; echo "$FILE"

        # create the source file
        cat template/int_nuclear.in > "$FILE"

        # calculate the expected energy
        V=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES -p v | ../script/extract/nuclear.sh)
        V=$(echo "$V" | tr "\n" " " | tr -s " " | sed "s/^ // ; s/\s*$//g ; s/ /, /g")

        # replace values in the source file
        sed -i "s/BNAME/$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/VVAL/$V/g" "$FILE"
    done
done

# Hartree-Fock
for SYSTEM in "${SYSTEMS[@]}"; do
    for BASIS in "${BASES[@]}"; do
        # echo the source file name
        FILE="energy_${SYSTEM,,}_hf_$BASIS.cpp"; echo "$FILE"

        # create the source file
        cat template/energy_hf.in > "$FILE"

        # calculate the expected energy
        ENERGY=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES hf -d 3 5 -m 1000 -t 1e-8 | grep FINAL | awk '{print $4}')

        # replace values in the source file
        sed -i "s/BNAME/$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/EVALUE/$ENERGY/g" "$FILE"
    done
done

# Hartree-Fock analytical gradient
for SYSTEM in "${SYSTEMS[@]}"; do
    for BASIS in "${BASES[@]}"; do
        # echo the source file name
        FILE="grad_${SYSTEM,,}_hf_$BASIS.cpp"; echo "$FILE"

        # create the source file
        cat template/grad_hf.in > "$FILE"

        # calculate the expected energy
        GRAD=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES hf -d 3 5 -m 1000 -t 1e-8 -g | sed -n "/ANAL/,/GRAD/p" | tail -n +5 | head -n -2)
        GRAD=$(echo "$GRAD" | tr "\n" " " | tr -s " " | sed "s/^ // ; s/\s*$//g ; s/ /, /g")

        # replace values in the source file
        sed -i "s/BNAME/$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/GRADVAL/$GRAD/g" "$FILE"
    done
done

# MP2
for SYSTEM in "${SYSTEMS[@]}"; do
    for BASIS in "${BASES[@]}"; do
        # echo the source file name
        FILE="energy_${SYSTEM,,}_mp2_$BASIS.cpp"; echo "$FILE"

        # create the source file
        cat template/energy_mp2.in > "$FILE"

        # calculate the expected energy
        ENERGY=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES hf -d 3 5 -m 1000 -t 1e-8 mp2 | grep "FINAL MP2" | awk '{print $4}')

        # replace values in the source file
        sed -i "s/BNAME/$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/EVALUE/$ENERGY/g" "$FILE"
    done
done

for FILE in *.cpp; do mv "$FILE" "$(echo "$FILE" | sed 's/*/s/g ; s/-//g')" &> /dev/null; done
