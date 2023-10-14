#!/bin/bash

HFOPT="-i 1000 -t 1e-8"; CORES=2

read -d '' HF << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/hazel -b \"%s\" -c %d -f \${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -n 2 -s %d hf $HFOPT)
set_tests_properties(%s PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")
EOF

read -d '' HFAG << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/hazel -b \"%s\" -c %d -f \${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -n 2 -s %d hf $HFOPT -g 0 1e-5)
set_tests_properties(%s PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")
EOF

read -d '' MP2 << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/hazel -b \"%s\" -c %d -f \${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -n 2 -s %d hf $HFOPT mp2)
set_tests_properties(%s PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")
EOF

create() {
    SYSTEM="$1"; CHARGE="$2"; MULT="$3"; BASIS="${4,,}"; METHOD="$5"; BASISFILE=$(echo "$BASIS" | sed -e 's/+/p/g ; s/*/s/g')

    if [ "$METHOD" == "HF" ]; then
        EXPECT=$(timeout 1s ./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -n $CORES -s "$MULT" hf $HFOPT | grep "FINAL HARTREE")
        [[ "$EXPECT" ]] && printf "\n\n# test-energy: %s %d %d %s HF\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"; TNAME="${SYSTEM}_${CHARGE}-${MULT}_${BASISFILE}_hf_energy"
        [[ "$EXPECT" ]] && printf "$HF" "$TNAME" "$BASIS" "$CHARGE" "$SYSTEM" "$MULT" "$TNAME" "$(echo ${EXPECT::-6} | sed -e 's/+/\\\\+/g')"
    elif [ "$METHOD" == "HFAG" ]; then
        EXPECT=$(timeout 1s ./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -n $CORES -s "$MULT" hf $HFOPT -g 0 1e-5 | grep "NORM")
        [[ "$EXPECT" ]] && printf "\n\n# test-gradient: %s %d %d %s HF\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"; TNAME="${SYSTEM}_${CHARGE}-${MULT}_${BASISFILE}_hf_gradient"
        [[ "$EXPECT" ]] && printf "$HFAG" "$TNAME" "$BASIS" "$CHARGE" "$SYSTEM" "$MULT" "$TNAME" "$(echo ${EXPECT} | sed -e 's/+/\\\\+/g')"
    elif [ "$METHOD" == "MP2" ]; then
        EXPECT=$(timeout 1s ./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -n $CORES -s "$MULT" hf $HFOPT mp2 | grep "FINAL MP2")
        [[ "$EXPECT" ]] && printf "\n\n# test-energy: %s %d %d %s MP2\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"; TNAME="${SYSTEM}_${CHARGE}-${MULT}_${BASISFILE}_mp2_energy"
        [[ "$EXPECT" ]] && printf "$MP2" "$TNAME" "$BASIS" "$CHARGE" "$SYSTEM" "$MULT" "$TNAME" "$(echo ${EXPECT::-6} | sed -e 's/+/\\\\+/g')"
    fi
}

printf '# TESTS %s' $(printf '=%.0s' {1..142})

for SYSTEM in \
"H2" "HF" "HCl" "water" "formaldehyde" "methane";
do
for BASIS in \
"mini" "midi" "sto-2g" "sto-3g" "sto-4g" "sto-5g" "sto-6g" "3-21g" "4-31g" "6-21g" "6-31g" "6-31+g" "6-31++g" "6-31+g*" "6-31++g*" "6-31+g**" "6-31++g**" "6-311g" "6-311+g" "6-311++g" "6-311+g*" \
"6-311++g*" "6-311+g**" "6-311++g**" "cc-pvdz" "cc-pvtz";
do
    create "$SYSTEM" 0 1 "$BASIS" "HF"; create "$SYSTEM" 0 1 "$BASIS" "MP2"; create "$SYSTEM" 0 1 "$BASIS" "HFAG"
    create "$SYSTEM" 1 2 "$BASIS" "HF"; create "$SYSTEM" -1 2 "$BASIS" "HF"; create "$SYSTEM" 2 1 "$BASIS" "HF"
done
done

echo ""
