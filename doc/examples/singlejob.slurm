#!/bin/bash --login

## Computational Resources to request:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=general

## "Wall time" (time limit), Days-HH:MM:SS
#SBATCH --time=0-1:00:00

## Optional job name
#SBATCH --job-name=SingleAmpsciJob

## Replace with your group's account string:
#SBATCH --account=a_your_account_string

module purge
module load foss/2024a
module load gsl/2.8-gcc-13.3.0

## Set number of OMP threads:
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
## Set BLAS threads to match (for OpenBLAS or MKL):
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
## If using MKL instead of OpenBLAS, replace above with:
## export MKL_THREADING_LAYER=GNU
## export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

## Path to where ampsci was compiled
AMPSCI=$HOME/ampsci/ampsci

## Input and output file names:
input=inputfile.in
output=${input/".in"/".out"}

$AMPSCI $input |tee -a $output
