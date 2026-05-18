#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --job-name=ArrayAmpsciJob
#SBATCH --time=0-1:00:00
#SBATCH --partition=general
## Replace with your group's account string:
#SBATCH --account=a_your_account

## Set array range to match number of input files (0-indexed):
#SBATCH --array=0-24

module purge
module load foss/2024a
module load gsl/2.8-gcc-13.3.0

## Set number of OMP threads:
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## Path to where ampsci was compiled
AMPSCI=$HOME/ampsci/ampsci

## Build list of input files (all matching "*.in" in current directory)
input_file_list=()
for input_file in *.in; do
  input_file_list+=($input_file)
done

## Input and output file names - one for each array job:
input=${input_file_list[${SLURM_ARRAY_TASK_ID}]}
output=${input/".in"/".out"}

$AMPSCI $input |tee -a $output
