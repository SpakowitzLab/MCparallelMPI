#/bin/bash

echo "Compile"
gfortran -ffpe-trap='' -g -O5 SIMcode/getpara.f95 \
                      MCcode/MC_Hamiltonian.f95\
                      MCcode/MC_interp.f95\
                      SIMcode/inputparams.f90\
                      SIMcode/simMod.f90\
                      MCcode/MC_int.f95\
                      r_2_phi.f90 -o r_to_phi 
mv r_to_phi ..
cd ..
echo "Now run"
# now run the output
rm data/phi_from_r*
./r_to_phi

cd code
