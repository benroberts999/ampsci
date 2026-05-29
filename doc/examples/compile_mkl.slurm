#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=9G
#SBATCH --job-name=CompileAMPSCI_MKL
#SBATCH --time=0:05:00
#SBATCH --partition=general
## Replace with your group's account string:
#SBATCH --account=a_your_account

## Uncomment and set for specific CPU architecture (Bunya etc.):
##SBATCH --constraint=epyc3

module purge
module load intel/2024a          # provides MKL
module load foss/2024a     # gets GCC + OpenBLAS (for GSL dependency)
module load gsl/2.8-gcc-13.3.0  # GSL (C library; compatible with GCC or Intel compiler)

./configure.sh -y

## -lmkl_rt :: link against MKL version of lapack/blas
make LDLIBS="-lgsl -lgslcblas -lmkl_rt" -j8 ampsci |tee -a compile-log-mkl.out
