#!/bin/bash

SYSTEMS=("ammonia" "ethane" "ethylene" "formaldehyde" "methane" "water")
BASES=("mini" "3-21g" "sto-3g" "6-31g" "6-31g*" "cc-pvdz")
INTEGRALS=("overlap" "kinetic" "nuclear")

HF="-d 3 5 -m 1000 -t 1e-8"

CORES=64

# integral matrices
for INTEGRAL in "${INTEGRALS[@]}"; do
    for SYSTEM in "${SYSTEMS[@]}"; do
        for BASIS in "${BASES[@]}"; do
            # echo the source file name
            FILE="integral_${INTEGRAL}_${SYSTEM,,}_$BASIS.cpp"; echo "$FILE"

            # create the source file
            cat "template/integral.in" > "$FILE"

            # calculate the expected energy
            I=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES ints -p all | "../script/extract/$INTEGRAL.sh")
            I=$(echo "$I" | tr "\n" " " | tr -s " " | sed "s/^ // ; s/\s*$//g ; s/ /, /g")

            # replace values in the source file
            sed -i "s/(int/_$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')(int/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
            sed -i "s/Integral::/Integral::${INTEGRAL^}/g" "$FILE" && sed -i "s/_INTEGRAL/_$INTEGRAL/g" "$FILE"
            sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/IVAL/$I/g" "$FILE"
        done
    done
done

# Hartree-Fock
for SYSTEM in "${SYSTEMS[@]}"; do
    for BASIS in "${BASES[@]}"; do
        # echo the source file name
        FILE="energy_hf_${SYSTEM,,}_$BASIS.cpp"; echo "$FILE"

        # create the source file
        cat template/energy_hf.in > "$FILE"

        # calculate the expected energy
        ENERGY=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES hf $HF | grep FINAL | awk '{print $4}')

        # replace values in the source file
        sed -i "s/(int/_$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')(int/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/EVALUE/$ENERGY/g" "$FILE"
    done
done

# Hartree-Fock analytical gradient
for SYSTEM in "${SYSTEMS[@]}"; do
    for BASIS in "${BASES[@]}"; do
        # echo the source file name
        FILE="gradient_hf_${SYSTEM,,}_$BASIS.cpp"; echo "$FILE"

        # create the source file
        cat template/gradient_hf.in > "$FILE"

        # calculate the expected energy
        GRAD=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES hf $HF -g | sed -n "/ANAL/,/GRAD/p" | tail -n +5 | head -n -2)
        GRAD=$(echo "$GRAD" | tr "\n" " " | tr -s " " | sed "s/^ // ; s/\s*$//g ; s/ /, /g")

        # replace values in the source file
        sed -i "s/(int/_$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')(int/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/GRADVAL/$GRAD/g" "$FILE"
    done
done

# MP2
for SYSTEM in "${SYSTEMS[@]}"; do
    for BASIS in "${BASES[@]}"; do
        # echo the source file name
        FILE="energy_mp2_${SYSTEM,,}_$BASIS.cpp"; echo "$FILE"

        # create the source file
        cat template/energy_mp2.in > "$FILE"

        # calculate the expected energy
        ENERGY=$(../bin/hazel -b "$BASIS" -f "../example/molecule/$SYSTEM.xyz" -n $CORES hf $HF mp2 | grep "FINAL MP2" | awk '{print $4}')

        # replace values in the source file
        sed -i "s/(int/_$(echo "$BASIS" | sed -e 's/*/s/g' -e 's/-//g')(int/g" "$FILE" && sed -i "s/BASIS/${BASIS^^}/g" "$FILE"
        sed -i "s/SYSTEM/$SYSTEM/g" "$FILE" && sed -i "s/EVALUE/$ENERGY/g" "$FILE"
    done
done

for FILE in *.cpp; do mv "$FILE" "$(echo "$FILE" | sed 's/*/s/g ; s/-//g')" &> /dev/null; done
