#/bin/bash
#sleep 30

echo "Compile"
rm MCparrll_out
cd code
# compile with mpi's fortran compiler
mpifort -c DATAcode/* mersenne_twister.f90
mpifort -c -fbounds-check -Wall -W -fmax-errors=5 -O5 SIMcode/*  MCcode/* 
mpifort *.o -o MCparrll_out 
rm *.o
cd ..
mv code/MCparrll_out .
mkdir -p data
mv data/* trash/

sleep 2
echo "Now run"
# now run the output
# --prefix used to avoid changing path
mpirun --prefix ~/openmpi/ -np 8 MCparrll_out


# cd code
# gfortran -O3 -fbounds-check -o wlcsim SIMcode/* BDcode/* DATAcode/* MCcode/*
# cd ..
# mv code/wlcsim .
# mv data/* trash/.
# ./wlcsim
