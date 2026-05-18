\page hpc HPC / SLURM
\ingroup getting_started
\brief Tips for running ampsci on HPC systems using SLURM

## HPC / SLURM

ampsci runs well on HPC systems via SLURM. It is a shared-memory (OpenMP)
code -- it does not use MPI -- so all parallelism is within a single node.

Useful references:

* [SLURM documentation](https://slurm.schedmd.com/documentation.html)
* [Bunya user guide](https://github.com/UQ-RCC/hpc-docs/blob/main/guides/Bunya-User-Guide.md)
* [Friday user guide](https://research.smp.uq.edu.au/friday-cluster/)

--------------------------------------------------------------------------------

## Key SLURM directives

| Directive | Meaning |
|-----------|---------|
| `--nodes=1` | single node (ampsci does not use MPI) |
| `--ntasks-per-node=1` | one task per node |
| `--cpus-per-task=N` | number of OpenMP threads; set to match `make -jN` and `OMP_NUM_THREADS` |
| `--mem=XG` | memory per node |
| `--time=D-H:MM:SS` | wall time (time limit) |
| `--partition=...` | queue/partition name (site-specific) |
| `--account=...` | account to charge (site-specific) |
| `--constraint=...` | request specific CPU architecture (optional, site-specific) |

`#!/bin/bash --login` is recommended -- it sources the user's login environment,
which ensures `module` is available.

--------------------------------------------------------------------------------

## Modules

Load required modules before running ampsci or `configure.sh`.
`module purge` first avoids conflicts from previously loaded modules.

<div class="shell-block">
```shell
module purge
module load foss/2024a gsl/2.8-gcc-13.3.0
module list
```
</div>

See [Compilation Details](\ref compilation) for more on module names and versions.

--------------------------------------------------------------------------------

## Example scripts

Three example SLURM scripts are provided in `doc/examples/`:

* \ref compile.slurm -- compile ampsci
* \ref singlejob.slurm -- run a single ampsci job
* \ref arrayjob.slurm -- run an array of jobs (e.g. over a parameter range)

### Compile job

\include compile.slurm

### Single job

\include singlejob.slurm

### Array job

\include arrayjob.slurm

--------------------------------------------------------------------------------

## Tips

* **Threads:** set `OMP_NUM_THREADS` to match `--cpus-per-task`, or ampsci will
  default to using all available cores on the node:

<div class="shell-block">
```shell
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./ampsci input.in
```
</div>

### Memory

* ampsci memory use depends on the basis size.
* 8--16 GB is often sufficient for small calculations; large MBPT calculations may need more.
* Use `ampsci -z <Basis>` to estimate memory requirements (can be very rough). e.g.,

<div class="shell-block">
```shell
  ./ampsci -z 35spdfgh
```
</div>

### Output

* redirect output to a file for later inspection. Useing `tee` is recommended
* `tee` will output to screen and to text file; `-a` means append:

<div class="shell-block">
```shell
./ampsci input.in |tee -a output.out
```
</div>
