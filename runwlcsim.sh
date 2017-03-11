#/bin/bash

echo "Compile"
rm MCparrll_out
cd code
# compile with mpi's fortran compiler
mpifort -c DATAcode/* mersenne_twister.f90
mpifort -c -fbounds-check -Wall -W -fmax-errors=5 -O1 SIMcode/*  MCcode/* 
mpifort *.o -o MCparrll_out 
rm *.o
cd ..
mv code/MCparrll_out .
mkdir -p data
mv data/* trash/

touch data/error
echo "Now run"
# now run the output
# --prefix used to avoid changing path
mpirun --oversubscribe -np 15 MCparrll_out
