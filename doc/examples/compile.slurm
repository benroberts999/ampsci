#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=9G
#SBATCH --job-name=CompileAMPSCI
#SBATCH --time=0:05:00
#SBATCH --partition=general
## Replace with your group's account string:
#SBATCH --account=a_your_account

## Uncomment and set for specific CPU architecture (Bunya etc.):
##SBATCH --constraint=epyc3

module purge
module load foss/2024a
module load gsl/2.8-gcc-13.3.0

./configure.sh -y

make -j8 ampsci |tee -a compile-log.out
