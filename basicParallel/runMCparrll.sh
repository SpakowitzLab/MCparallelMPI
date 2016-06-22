#!/bin/sh
# compile with mpi's fortran compiler
mpifort mersenne_twister.f90 MCparrll_mpi.f90 -o MCparrll_out 

# remove all old data
rm data/*

# now run the output
# --prefix used to avoid changing path
mpirun --prefix ~/openmpi/ -np 4 MCparrll_out
