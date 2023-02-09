#!/bin/bash

BASIS="STO-3G"
METHOD="RHF"

[[ -f clean.sh ]] && ./clean.sh

echo "#!/bin/bash
rm -f clean.sh *.0 *.densities *.engrad *.gbw *.ges *.hess *.inp *.opt *.out *.tmp *.txt *.xyz *trj* " > clean.sh
chmod +x clean.sh

cat > water.xyz <<- EOM
3
water
O -0.066989 -0.023205 -0.026795
H  0.684834 -0.612063  0.519104
H  0.378371  0.980368 -0.093822
EOM

for MOL in *.xyz; do
    echo -e "! $METHOD $BASIS OPT LARGEPRINT\n*xyzfile 0 1 $MOL" > "${MOL%.*}.inp" && orca "${MOL%.*}.inp" | tee "${MOL%.*}.out"
    ENERGY=$(grep "FINAL SINGLE POINT ENERGY" "${MOL%.*}.out" | tail -n 1 | awk '{print $5}')
    sed -i "2s/.*/${MOL%.*}\/$METHOD\/$BASIS\/$ENERGY/" "$MOL"
done
