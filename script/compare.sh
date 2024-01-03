echo "System,Basis,NBF,Method,Hazel Energy,Hazel Time,Orca Energy,Orca Time,Orca Energy Difference"

# define the systems to compare
SYSTEMS=("HCl" "CO2" "water" "ammonia" "formaldehyde" "methane" "ethylene" "cyclopropanone" "ethanol" "acetone" "pyrrole" "benzene")

# define the bases to use
BASES=("mini" "sto-3g" "6-311++g**" "cc-pvdz" "def2-svp" "aug-cc-pvdz")

# compare RHF
for SYSTEM in "${SYSTEMS[@]}"; do
for BASIS in "${BASES[@]}"; do
    NBF=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" | grep NBF | awk '{print $7}'); [[ $NBF -ge 80 ]] && continue
    echo -n "$SYSTEM,$BASIS,$NBF,HF,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8)
    HAZELE=$(echo "$HAZEL" | grep "FINAL HA" | awk '{print $4}' | head -c -3)
    HAZELT=$(echo "$HAZEL" | grep "TOTAL EX" | awk '{print $4}')
    echo -n "$HAZELE,$HAZELT,"
    ORCA=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" orca -m rhf)
    ORCAE=$(echo "$ORCA" | grep "FINAL EN" | awk '{print $3}' | head -c -3)
    ORCAT=$(echo "$ORCA" | grep "TOTAL EX" | awk '{print $4}')
    echo -n "$ORCAE,$ORCAT,"
    DIFF=$(echo "($HAZELE)-($ORCAE)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done

# compare MP2
for SYSTEM in "${SYSTEMS[@]}"; do
for BASIS in "${BASES[@]}"; do
    NBF=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" | grep NBF | awk '{print $7}'); [[ $NBF -ge 50 ]] && continue
    echo -n "$SYSTEM,$BASIS,$NBF,MP2,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 mp2)
    HAZELE=$(echo "$HAZEL" | grep "FINAL MP" | awk '{print $4}' | head -c -3)
    HAZELT=$(echo "$HAZEL" | grep "TOTAL EX" | awk '{print $4}')
    echo -n "$HAZELE,$HAZELT,"
    ORCA=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" orca -m mp2)
    ORCAE=$(echo "$ORCA" | grep "FINAL EN" | awk '{print $3}' | head -c -3)
    ORCAT=$(echo "$ORCA" | grep "TOTAL EX" | awk '{print $4}')
    echo -n "$ORCAE,$ORCAT,"
    DIFF=$(echo "($HAZELE)-($ORCAE)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done

# compare CISD
for SYSTEM in "${SYSTEMS[@]}"; do
for BASIS in "${BASES[@]}"; do
    NBF=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" | grep NBF | awk '{print $7}'); [[ $NBF -ge 10 ]] && continue
    echo -n "$SYSTEM,$BASIS,$NBF,CISD,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 cisd)
    HAZELE=$(echo "$HAZEL" | grep "FINAL CI" | awk '{print $4}' | head -c -3)
    HAZELT=$(echo "$HAZEL" | grep "TOTAL EX" | awk '{print $4}')
    echo -n "$HAZELE,$HAZELT,"
    ORCA=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" orca -m cisd)
    ORCAE=$(echo "$ORCA" | grep "FINAL EN" | awk '{print $3}' | head -c -3)
    ORCAT=$(echo "$ORCA" | grep "TOTAL EX" | awk '{print $4}')
    echo -n "$ORCAE,$ORCAT,"
    DIFF=$(echo "($HAZELE)-($ORCAE)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done

# compare FCI
for SYSTEM in "${SYSTEMS[@]}"; do
for BASIS in "${BASES[@]}"; do
    NBF=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" | grep NBF | awk '{print $7}'); [[ $NBF -ge 10 ]] && continue
    echo -n "$SYSTEM,$BASIS,$NBF,FCI,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 fci)
    HAZELE=$(echo "$HAZEL" | grep "FINAL CI" | awk '{print $4}' | head -c -3)
    HAZELT=$(echo "$HAZEL" | grep "TOTAL EX" | awk '{print $4}')
    echo -n "$HAZELE,$HAZELT,"
    ORCA=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" orca -m fci)
    ORCAE=$(echo "$ORCA" | grep "FINAL EN" | awk '{print $3}' | head -c -3)
    ORCAT=$(echo "$ORCA" | grep "TOTAL EX" | awk '{print $4}')
    echo -n "$ORCAE,$ORCAT,"
    DIFF=$(echo "($HAZELE)-($ORCAE)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done
