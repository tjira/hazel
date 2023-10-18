#!/bin/bash

HFOPT="-i 1000 -t 1e-8"; CORES=2

read -d '' RHF << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/hazel -b \"%s\" -c %d -f \${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -n 2 -s %d rhf $HFOPT)
set_tests_properties(%s PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")
EOF

read -d '' RHFAG << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/hazel -b \"%s\" -c %d -f \${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -n 2 -s %d rhf $HFOPT -g 0 1e-5)
set_tests_properties(%s PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")
EOF

read -d '' RMP2 << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/hazel -b \"%s\" -c %d -f \${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -n 2 -s %d rhf $HFOPT mp2)
set_tests_properties(%s PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")
EOF

read -d '' UHF << EOF
add_test(NAME %s COMMAND \${PROJECT_SOURCE_DIR}/bin/hazel -b \"%s\" -c %d -f \${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -n 2 -s %d uhf $HFOPT)
set_tests_properties(%s PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")
EOF

create() {
    SYSTEM="$1"; CHARGE="$2"; MULT="$3"; BASIS="${4,,}"; METHOD="$5"; BASISFILE=$(echo "$BASIS" | sed -e 's/+/p/g ; s/*/s/g')

    if [ "$METHOD" == "RHF" ]; then
        EXPECT=$(./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -n $CORES -s "$MULT" rhf $HFOPT | grep "FINAL HARTREE")
        [[ "$EXPECT" ]] && printf "\n\n# test-energy: %s %d %d %s RHF\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"; TNAME="${SYSTEM}_${CHARGE}-${MULT}_${BASISFILE}_rhf_energy"
        [[ "$EXPECT" ]] && printf "$RHF" "$TNAME" "$BASIS" "$CHARGE" "$SYSTEM" "$MULT" "$TNAME" "$(echo ${EXPECT::-6} | sed -e 's/+/\\\\+/g')"
    elif [ "$METHOD" == "RHFAG" ]; then
        EXPECT=$(./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -n $CORES -s "$MULT" rhf $HFOPT -g 0 1e-5 | grep "NORM")
        [[ "$EXPECT" ]] && printf "\n\n# test-gradient: %s %d %d %s RHF\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"; TNAME="${SYSTEM}_${CHARGE}-${MULT}_${BASISFILE}_rhf_gradient"
        [[ "$EXPECT" ]] && printf "$RHFAG" "$TNAME" "$BASIS" "$CHARGE" "$SYSTEM" "$MULT" "$TNAME" "$(echo ${EXPECT} | sed -e 's/+/\\\\+/g')"
    elif [ "$METHOD" == "RMP2" ]; then
        EXPECT=$(./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -n $CORES -s "$MULT" rhf $HFOPT mp2 | grep "FINAL MP2")
        [[ "$EXPECT" ]] && printf "\n\n# test-energy: %s %d %d %s RMP2\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"; TNAME="${SYSTEM}_${CHARGE}-${MULT}_${BASISFILE}_rmp2_energy"
        [[ "$EXPECT" ]] && printf "$RMP2" "$TNAME" "$BASIS" "$CHARGE" "$SYSTEM" "$MULT" "$TNAME" "$(echo ${EXPECT::-6} | sed -e 's/+/\\\\+/g')"
    elif [ "$METHOD" == "UHF" ]; then
        EXPECT=$(./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -n $CORES -s "$MULT" uhf $HFOPT | grep "FINAL HARTREE")
        [[ "$EXPECT" ]] && printf "\n\n# test-energy: %s %d %d %s UHF\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"; TNAME="${SYSTEM}_${CHARGE}-${MULT}_${BASISFILE}_rhf_energy"
        [[ "$EXPECT" ]] && printf "$UHF" "$TNAME" "$BASIS" "$CHARGE" "$SYSTEM" "$MULT" "$TNAME" "$(echo ${EXPECT::-6} | sed -e 's/+/\\\\+/g')"
    fi
}

printf '# TESTS %s' $(printf '=%.0s' {1..142})

for SYSTEM in \
"water";
do
for BASIS in \
"mini" "midi" "sto-2g" "sto-3g" "sto-4g" "sto-5g" "sto-6g" "3-21g" "4-31g" "6-21g" "6-31g" "6-31+g" "6-31++g" "6-31+g*" "6-31++g*" "6-31+g**" "6-31++g**" "6-311g" "6-311+g" "6-311++g" "6-311+g*" \
"6-311++g*" "6-311+g**" "6-311++g**" "cc-pvdz" "cc-pvtz";
do
    create "$SYSTEM" 0 1 "$BASIS" "RHF"; create "$SYSTEM" 0 1 "$BASIS" "RMP2"; create "$SYSTEM" 0 1 "$BASIS" "RHFAG"
    create "$SYSTEM" 1 2 "$BASIS" "UHF"; create "$SYSTEM" -1 2 "$BASIS" "UHF"; create "$SYSTEM" 2 1 "$BASIS" "UHF"
done
done

echo ""
