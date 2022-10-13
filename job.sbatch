#!/bin/bash
#SBATCH -J rbtest
#SBATCH -o test.o%j
#SBATCH -t 48:00:00
#SBATCH -N 2 -n 52
#SBATCH -A ostilla
#SBATCH --mail-type=END
#SBATCH --mail-user=rostilla@central.uh.edu
#SBATCH -p batch

module load HDF5/1.12.1-iimpi-2022.00 FFTW/3.3.10-iimpi-2022.00 intel/.2022.00

ulimit -s unlimited
mpirun  -np 52 ./boutnp
