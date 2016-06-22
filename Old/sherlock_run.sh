#!/bin/bash
#SBATCH --job-name=L20N10_a2
#SBATCH --output=SOut.L20N10_a2
#SBATCH --error=SError.L20N10_a2
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --mem=2000
#SBATCH --ntasks-per-node=16
./wlcsim pt-1
