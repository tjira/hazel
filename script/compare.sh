echo "system,basis,nbf,method,hazel,orca,difference"

# define the systems to compare
SYSTEMS=("HCl" "CO2" "water" "ammonia" "formaldehyde" "methane" "ethylene" "cyclopropanone" "ethanol" "acetone" "pyrrole" "benzene")

# define the bases to use
BASES=("mini" "sto-3g" "6-311++g**" "cc-pvdz" "def2-svp" "aug-cc-pvdz")

# compare RHF
for SYSTEM in "${SYSTEMS[@]}"; do
for BASIS in "${BASES[@]}"; do
    NBF=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" | grep NBF | awk '{print $7}'); [[ $NBF -ge 80 ]] && continue
    echo -n "$SYSTEM,$BASIS,$NBF,HF,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 | grep "FINAL" | awk '{print $4}' | head -c -3)
    echo -n "$HAZEL,"
    ORCA=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" orca -m rhf | grep "FINAL" | awk '{print $3}' | head -c -3)
    echo -n "$ORCA,"
    DIFF=$(echo "($HAZEL)-($ORCA)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done

# compare MP2
for SYSTEM in "${SYSTEMS[@]}"; do
for BASIS in "${BASES[@]}"; do
    NBF=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" | grep NBF | awk '{print $7}'); [[ $NBF -ge 50 ]] && continue
    echo -n "$SYSTEM,$BASIS,$NBF,MP2,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 mp2 | grep "FINAL MP2" | awk '{print $4}' | head -c -3)
    echo -n "$HAZEL,"
    ORCA=$(./script/orca.sh "example/molecule/$SYSTEM.xyz" 0 1 mp2 "$BASIS" | head -n 1 | awk '{print $2}')
    echo -n "$ORCA,"
    DIFF=$(echo "($HAZEL)-($ORCA)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done

# compare CISD
for SYSTEM in "${SYSTEMS[@]}"; do
for BASIS in "${BASES[@]}"; do
    NBF=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" | grep NBF | awk '{print $7}'); [[ $NBF -ge 10 ]] && continue
    echo -n "$SYSTEM,$BASIS,$NBF,CISD,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 cisd | grep "FINAL CI" | awk '{print $4}' | head -c -3)
    echo -n "$HAZEL,"
    ORCA=$(./script/orca.sh "example/molecule/$SYSTEM.xyz" 0 1 cisd "$BASIS" | head -n 1 | awk '{print $2}')
    echo -n "$ORCA,"
    DIFF=$(echo "($HAZEL)-($ORCA)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done

# compare FCI
for SYSTEM in "${SYSTEMS[@]}"; do
for BASIS in "${BASES[@]}"; do
    NBF=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" | grep NBF | awk '{print $7}'); [[ $NBF -ge 10 ]] && continue
    echo -n "$SYSTEM,$BASIS,$NBF,FCI,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 fci | grep "FINAL CI" | awk '{print $4}' | head -c -3)
    echo -n "$HAZEL,"
    ORCA=$(./script/orca.sh "example/molecule/$SYSTEM.xyz" 0 1 fci "$BASIS" | head -n 1 | awk '{print $2}')
    echo -n "$ORCA,"
    DIFF=$(echo "($HAZEL)-($ORCA)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done
