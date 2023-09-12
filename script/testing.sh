#!/bin/bash

SYSTEMS=("ammonia" "ethylene" "methane" "water")
BASES=("3-21g" "sto-3g" "6-31g" "6-31g*")
CMS=("0:1" "1:2" "-1:2")

printf '# TESTS %s\n' $(printf '=%.0s' {1..142})

# Hartree-Fock Method
for SYSTEM in "${SYSTEMS[@]}"; do
    for CM in "${CMS[@]}"; do
        for BASIS in "${BASES[@]}"; do
            # extract charge and multiplicity
            CMSPLIT=(${CM//:/ }); CHARGE=${CMSPLIT[0]}; MULT=${CMSPLIT[1]}

            # print the test specification
            printf "\n# test-energy: %s %d %d %s HF\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"

            # extract the expected output
            EXPECT=$(./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -s "$MULT" hf -i 1000 -t 1e-8 | grep "FINAL HARTREE")

            # print the test commands
            printf 'add_test(NAME %s_%d-%d_%s_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "%s" -c %d -f ${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -s %d hf -i 1000 -t 1e-8)\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "$BASIS" "$CHARGE" "$SYSTEM" "$MULT"
            printf 'set_tests_properties(%s_%d-%d_%s_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "${EXPECT::-6}"
        done
    done
done

# Møller–Plesset of Second Order
for SYSTEM in "${SYSTEMS[@]}"; do
    for CM in "${CMS[@]:0:1}"; do
        for BASIS in "${BASES[@]}"; do
            # extract charge and multiplicity
            CMSPLIT=(${CM//:/ }); CHARGE=${CMSPLIT[0]}; MULT=${CMSPLIT[1]}

            # print the test specification
            printf "\n# test-energy: %s %d %d %s MP2\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"

            # extract the expected output
            EXPECT=$(./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -s "$MULT" hf -i 1000 -t 1e-8 mp2 | grep "FINAL MP2")

            # print the test commands
            printf 'add_test(NAME %s_%d-%d_%s_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "%s" -c %d -f ${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -s %d hf -i 1000 -t 1e-8 mp2)\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "$BASIS" "$CHARGE" "$SYSTEM" "$MULT"
            printf 'set_tests_properties(%s_%d-%d_%s_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "${EXPECT::-6}"
        done
    done
done

# Hartree-Fock Analytical Gradient
for SYSTEM in "${SYSTEMS[@]}"; do
    for CM in "${CMS[@]:0:1}"; do
        for BASIS in "${BASES[@]}"; do
            # extract charge and multiplicity
            CMSPLIT=(${CM//:/ }); CHARGE=${CMSPLIT[0]}; MULT=${CMSPLIT[1]}

            # print the test specification
            printf "\n# test-angrad: %s %d %d %s HF\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"

            # extract the expected output
            EXPECT=$(./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -s "$MULT" hf -i 1000 -t 1e-8 -g 0 1e-5 | grep "GRADIENT NORM")

            # print the test commands
            printf 'add_test(NAME %s_%d-%d_%s_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "%s" -c %d -f ${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -s %d hf -i 1000 -t 1e-8 -g 0 1e-5)\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "$BASIS" "$CHARGE" "$SYSTEM" "$MULT"
            printf 'set_tests_properties(%s_%d-%d_%s_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "${EXPECT}"
        done
    done
done

# Hartree-Fock Numerical Gradient
for SYSTEM in "${SYSTEMS[@]}"; do
    for CM in "${CMS[@]}"; do
        for BASIS in "${BASES[@]}"; do
            # extract charge and multiplicity
            CMSPLIT=(${CM//:/ }); CHARGE=${CMSPLIT[0]}; MULT=${CMSPLIT[1]}

            # print the test specification
            printf "\n# test-numgrad: %s %d %d %s HF\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"

            # extract the expected output
            EXPECT=$(./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -s "$MULT" hf -i 1000 -t 1e-8 -g 1 1e-5 | grep "GRADIENT NORM")

            # print the test commands
            printf 'add_test(NAME %s_%d-%d_%s_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "%s" -c %d -f ${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -s %d hf -i 1000 -t 1e-8 -g 1 1e-5)\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "$BASIS" "$CHARGE" "$SYSTEM" "$MULT"
            printf 'set_tests_properties(%s_%d-%d_%s_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "${EXPECT}"
        done
    done
done

# Møller–Plesset of Second Order Numerical Gradient
for SYSTEM in "${SYSTEMS[@]}"; do
    for CM in "${CMS[@]:0:1}"; do
        for BASIS in "${BASES[@]}"; do
            # extract charge and multiplicity
            CMSPLIT=(${CM//:/ }); CHARGE=${CMSPLIT[0]}; MULT=${CMSPLIT[1]}

            # print the test specification
            printf "\n# test-numgrad: %s %d %d %s MP2\n" "$SYSTEM" "$CHARGE" "$MULT" "${BASIS^^}"

            # extract the expected output
            EXPECT=$(./bin/hazel -b "$BASIS" -c "$CHARGE" -f "./example/molecule/$SYSTEM.xyz" -s "$MULT" hf -i 1000 -t 1e-8 mp2 -g 1 1e-5 | grep "GRADIENT NORM")

            # print the test commands
            printf 'add_test(NAME %s_%d-%d_%s_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "%s" -c %d -f ${PROJECT_SOURCE_DIR}/example/molecule/%s.xyz -s %d hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "$BASIS" "$CHARGE" "$SYSTEM" "$MULT"
            printf 'set_tests_properties(%s_%d-%d_%s_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "%s")\n' \
                   "$SYSTEM" "$CHARGE" "$MULT" $(echo "$BASIS" | sed -e 's/+/p/g' -e 's/*/s/g') "${EXPECT}"
        done
    done
done
