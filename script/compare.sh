echo "system,basis,method,hazel,orca,difference"

for SYSTEM in "water" "formaldehyde" "acetone"; do
for BASIS in "mini" "midi" "sto-3g" "3-21g" "6-31g" "6-31+g*" "6-31+g**" "6-31++g**" "6-311g" "6-311+g*" "6-311+g**" "6-311++g**" "cc-pvdz"; do
    echo -n "$SYSTEM,$BASIS,HF,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 | grep "FINAL" | awk '{print $4}' | head -c -3)
    echo -n "$HAZEL,"
    ORCA=$(./script/orca.sh "example/molecule/$SYSTEM.xyz" 0 1 hf "$BASIS" | head -n 1 | awk '{print $2}')
    echo -n "$ORCA,"
    DIFF=$(echo "($HAZEL)-($ORCA)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done

for SYSTEM in "water" "formaldehyde" "acetone"; do
for BASIS in "mini" "midi" "sto-3g" "3-21g" "6-31g" "6-31+g*" "6-31+g**" "6-31++g**" "6-311g" "6-311+g*" "6-311+g**" "6-311++g**" "cc-pvdz"; do
    echo -n "$SYSTEM,$BASIS,MP2,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 mp2 | grep "FINAL MP2" | awk '{print $4}' | head -c -3)
    echo -n "$HAZEL,"
    ORCA=$(./script/orca.sh "example/molecule/$SYSTEM.xyz" 0 1 mp2 "$BASIS" | head -n 1 | awk '{print $2}')
    echo -n "$ORCA,"
    DIFF=$(echo "($HAZEL)-($ORCA)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done

for SYSTEM in "H2" "HF" "water"; do
for BASIS in "mini" "sto-3g"; do
    echo -n "$SYSTEM,$BASIS,FCI,"
    HAZEL=$(./bin/hazel -b "$BASIS" -f "example/molecule/$SYSTEM.xyz" rhf -i 1000 -t 1e-8 fci | grep "FINAL CI" | awk '{print $4}' | head -c -3)
    echo -n "$HAZEL,"
    ORCA=$(./script/orca.sh "example/molecule/$SYSTEM.xyz" 0 1 fci "$BASIS" | head -n 1 | awk '{print $2}')
    echo -n "$ORCA,"
    DIFF=$(echo "($HAZEL)-($ORCA)" | bc -l)
    printf "%.2e\n" "${DIFF#-}"
done
done
