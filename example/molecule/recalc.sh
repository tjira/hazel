#!/bin/bash

BASIS="STO-3G"
METHOD="RHF"

echo "#!/bin/bash
rm -f clean.sh *.0 *.densities *.engrad *.gbw *.hess *.inp *.opt *.out *.tmp *.txt *trj* " > clean.sh
chmod +x clean.sh

for MOL in *.xyz; do
    echo -e "! $METHOD $BASIS LARGEPRINT\n*xyzfile 0 1 $MOL" > "${MOL%.*}.inp" && orca "${MOL%.*}.inp" | tee "${MOL%.*}.out"
done
