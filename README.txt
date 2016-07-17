#-------------------------------------------------------------------------------------
#
#     Basic info
#
#-------------------------------------------------------------------------------------
Authors: 
Andrew Spakowitz  (wrote initial code for wlc)
Sifan Mao (Improved MC setup)
Quinn MacPherson (Improved MC setup and MPI paralleization)  spring 2016
Contact: qmac@stanford.edu

Code diagram:
https://drive.google.com/file/d/0B9VKZRfyscD5WGhqWkFabUhmWWs/view?usp=sharing


#---------------------------------------------------------------------------------------
#
#   Example running script
#
#------------------------------------------------------------------------------------------

#/bin/bash
sleep 30

echo "Compile"

cd code
# compile with mpi's fortran compiler
mpifort -fbounds-check -O3 mersenne_twister.o SIMcode/* DATAcode/* MCcode/* -o MCparrll_out 
cd ..
mv code/MCparrll_out .
mv data/* trash/


echo "Now run"
# now run the output
# --prefix used to avoid changing path
mpirun --prefix ~/openmpi/ -np 40 MCparrll_out


# cd code
# gfortran -O3 -fbounds-check -o wlcsim SIMcode/* BDcode/* DATAcode/* MCcode/*
# cd ..
# mv code/wlcsim .
# mv data/* trash/.
# ./wlcsim


# -----------------------------------------------------------------------------------
#
#   How to make job non-interactive
#
# ----------------------------------------------------------------------------------

# Start the job interactively
./runwlcsim  # example script

# While sleeping and before it splits to many threads
[ctrl] z
disown -h %1  # or a different number
bg 1
logout


jobs # show jobs
top  # show computer use

# ------------------------------------------------------------------------------------
#
#   Installing MPI on tower
#
# ------------------------------------------------------------------------------------
#Downloaded mpi to laptop
#Tranfered to my home dir

# install

tar -zxvf openmpi...tar.gz
./configure --prefix=/usr/local
sudo make all install

# add to path

vi ~/.bash_profile
PATH=$PATH:/usr/local/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
alis ls='ls --color'
