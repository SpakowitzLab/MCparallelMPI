#/bin/bash

rm snap0* 2>/dev/null
rm temp 2>/dev/null
touch temp
for i in `seq 15 34`; 
    do name=$(printf "snap%04d.pdb" $i); 
    head -1000 ../../K7812/data/r$[i]v2 | tail -20 >> temp;
done

python r2pdbTime.py temp > time.pdb;

pymol mkTime.pml
