#/bin/bash

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

echo "Now run"
# now run the output
# --prefix used to avoid changing path
mpirun -np 8 MCparrll_out


