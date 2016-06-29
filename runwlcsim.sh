#/bin/bash
#sleep 30

echo "Compile"

cd code
# compile with mpi's fortran compiler
mpifort -fbounds-check -O3 mersenne_twister.o SIMcode/* DATAcode/* MCcode/* -o MCparrll_out 
cd ..
mv code/MCparrll_out .
mv data/* trash/

sleep 2
echo "Now run"
# now run the output
# --prefix used to avoid changing path
mpirun --prefix ~/openmpi/ -np 4 MCparrll_out


# cd code
# gfortran -O3 -fbounds-check -o wlcsim SIMcode/* BDcode/* DATAcode/* MCcode/*
# cd ..
# mv code/wlcsim .
# mv data/* trash/.
# ./wlcsim
