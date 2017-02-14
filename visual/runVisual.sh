#/bin/bash

rm snap0*

for i in `seq 1 14`; 
  do name=$(printf "snap%04d.pdb" $i); 
  python r2pdb.py ../../K312/data/r29v$i > "${name}";
done

pymol mksnap.pml
