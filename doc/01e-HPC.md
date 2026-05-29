\page hpc HPC / SLURM

\brief Tips for running ampsci on HPC systems using SLURM

## HPC / SLURM

ampsci runs well on HPC systems via SLURM.
It is a shared-memory (OpenMP)
code (it does not use MPI) so all parallelism is within a single node.

On HPC systems, compilation and jobs are typically submitted via SLURM -- see below for example scripts.
Do not run jobs, including compiling ampsci, directly on the login node of the HPC. You must submit all jobs via the queueing system.

* SLURM is a widely-used open-source job scheduler for HPC clusters; other schedulers exist (PBS, LSF, etc.) but SLURM is the most common.
* Jobs are submitted to a queue via `sbatch`, and then run when the required resources are available
  * The job script contains `#SBATCH` directives specifying resource requests, followed by the shell commands to run
  * You must request resources carefully: too little memory and your job will be killed; too much and you'll wait longer in the queue
  * Several example slurm job scripts are given below
* You can monitor the queue with `squeue`, and cancel jobs with `scancel`.

Useful references:

* [SLURM documentation](https://slurm.schedmd.com/documentation.html)
* [Bunya user guide](https://github.com/UQ-RCC/hpc-docs/blob/main/guides/Bunya-User-Guide.md)
* [Friday user guide](https://research.smp.uq.edu.au/friday-cluster/)

### Key SLURM commands

| Command | Description |
|---------|-------------|
| `sbatch job.slurm` | submit a job script to the queue |
| `squeue -u $USER` | list your queued/running jobs |
| `scancel <jobid>` | cancel a job |
| `sinfo` | show available partitions and node status |
| `sacct -j <jobid>` | show accounting info for a completed job |
| `seff <jobid>` | show CPU and memory efficiency for a completed job |

### Key SLURM directives

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

## Compiling ampsci on HPC systems

* On HPC systems, you will typically need to `module load` the required dependencies
  * What we require: C++ compiler, lapack, blas, GSL
  * Often, most of these come 'bundled' in a "toolchain" (e.g., foss)
* Load the required modules _before_ running `configure.sh`.
* On most systems (e.g., Friday) unversioned names work:

<div class="shell-block">
```shell
module load foss gsl
./configure.sh -y
make
```
</div>

On others (e.g. Bunya), explicit versions are required and must be consistent.
On most systems, the GSL module typically matches the foss toolchain something like (though specifics may change on different HPC systems):

| foss    | GCC    | Matching GSL (typical)|
|---------|--------|-----------------------|
| 2022a   | 11.3.0 | gsl/2.7-gcc-11.3.0    |
| 2023a   | 12.3.0 | gsl/2.7-gcc-12.3.0    |
| 2024a   | 13.3.0 | gsl/2.8-gcc-13.3.0    |

So, we would do something like:

<div class="shell-block">
```shell
module load foss/2024a gsl/2.8-gcc-13.3.0
./configure.sh -y
make
```
</div>

Load required modules before running ampsci or `configure.sh`.
It might be a good idea to add a `module purge` first, which avoids conflicts from previously loaded modules.

<div class="shell-block">
```shell
module purge
module load foss/2024a gsl/2.8-gcc-13.3.0
module list
```
</div>

Use `module avail foss` or `module spider gsl` or similar to find available versions.

`configure.sh` attempts to auto-detect BLAS/OpenBLAS version (via `$EBROOTOPENBLAS`) and sets `LDLIBS` accordingly.
If auto-detection fails, you will have to set it manually in `Makefile`, e.g.:

```makefile
LDLIBS ?= -lgsl -lgslcblas -lopenblas
```

Note: `-lgfortran` may also be required on some older configurations.

If `configure.sh` does not produce a working build, refer to the manual [Compilation Details](\ref compilation) for details

### BLAS threading (recommended for CI calculations)

The `foss` toolchain includes OpenBLAS, which is already multi-threaded.
Set `OPENBLAS_NUM_THREADS` in your job script to match `--cpus-per-task`:

```bash
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
```

**Intel MKL** is an alternative to OpenBLAS and may be faster on Intel nodes.
Load the MKL module (name is site-specific -- try `module spider mkl` or `module spider imkl`):

```bash
module load imkl   # name varies: intel-mkl, imkl, mkl, ...
```

Then set in the Makefile:

```makefile
LDLIBS ?= -lgsl -lgslcblas -lmkl_rt
```

And in your job script:

```bash
export MKL_THREADING_LAYER=GNU          # required with GNU OpenMP
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
```

See [Compilation Details](\ref compilation) for more on BLAS options.

--------------------------------------------------------------------------------

## Example scripts

Four example SLURM scripts are provided in `doc/examples/`:

* \ref compile.slurm -- compile ampsci (OpenBLAS)
* \ref compile_mkl.slurm -- compile ampsci with Intel MKL
* \ref singlejob.slurm -- run a single ampsci job
* \ref arrayjob.slurm -- run an array of jobs (e.g. over a parameter range)

### Compile job

\include compile.slurm

### Compile job (Intel MKL)

\include compile_mkl.slurm

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
