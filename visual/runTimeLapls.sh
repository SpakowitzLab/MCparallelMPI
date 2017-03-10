#/bin/bash

rm snap0* 2>/dev/null
rm temp 2>/dev/null
touch temp
for i in `seq 80 90`; 
    do name=$(printf "snap%04d.pdb" $i); 
    head -1000 ../Jie3p4sPEG1500n4HfPEG50/data/r$[i]v78 | tail -50 >> temp;
# Jie3p4sPEG1500n4HfPEG50
# ../Jie1p47sPEG900n2
# ../Jie2p28sPEG900n2
# ../Jie3p4sPEG1500n2
# ../Jie4p57sPEG1500n2
# ../Jie1p47sPEG900N
# ../Jie2p28sPEG900N
# ../Jie3p4sPEG1500N
# ../Jie4p57sPEG1500N
done

python r2pdbTime.py temp > time.pdb;

pymol mkTime.pml
