# Relativistic, self-consistent atomic structure code.

Solves the Dirac equation for atomic systems using the Hartree-Fock method.
Fully relativistic, includes finite-nuclear size, and can
solve for continuum states (energy normalisation).

 * With reasonable choices for the integration grids, typically converges
to better than a few parts in 10^16

 * Includes an option to vary the effective speed of light -
allowing non-relativistic approximation.

 * Wavefunctions are in form psi = (1/r) [f,ig], (using Dirac basis)

The part that solves the Dirac eigenvalue DE is based on book by W. Johnson,
with a few extensions that improve numerical stability and accuracy
 [W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007).]

### Compiling and use:

 * All programs compiled using the Makefile
 (run _$make_ or _$make programName_)
 * Must have LAPACK, GSL libraries installed already (see below)
 * Must create a directory called _./obj/_
   * code places object files inside here
 * All executables run like _$./programName_. By default, placed in same
  directory as makefile (can change this in make)
 * All programs have input options, stored and read from 'programName.in' file
   (can also give a different input file on runtime, e.g., to read input from
    `otherFile.txt`: _$./programName otherFile.txt_)
 * Note: below just tells how to use existing programs, to see how they work,
 see the comments/instructions inside the source code (all in /src/)
 * Tested with g++ and clang++.

Note: makes use of GSL libraries: https://www.gnu.org/software/gsl/

  * For example, with ubuntu: _$sudo apt-get install libgsl-dev_
  * Also needs LAPACK/BLAS libraries:
  _$sudo apt-get install libatlas-base-dev_ [and _liblapack-dev_ , but not _libblas-dev_ ?]

The above instructions are for linux (ubuntu). For windows, the easiest way (for me, anyway) is to make use of the recent 'windows subsystem for linux'. Instructions on installation/use here: https://www.roberts999.com/posts/2018/11/wsl-coding-windows-ubuntu
Then, the compilation + use can proceed as per above.

## hartreeFock

 * Solves relativistic Hartree Fock potential for core + valence states
 * Takes in core configuration: [Noble gas],extra (comma separated, no spaces)
 * (As well as Noble gas, can use Zn,Cd,Hg,Cn)
 * Can also add negative values: e.g.,

E.g. (V^N-1):
   * For Cs: '[Xe]'
   * For Au: '[Xe],4f14,5d10' OR '[Hg],6s-2'
   * For Tl: '[Xe],4f14,5d10,6s2' OR '[Hg]'
   * For I (V^N): '[Cd],5p5' OR '[Xe],5p-1'


## atomicKernal

 * Calculates the "Atomic Kernal" (for scattering/ionisation) for each core
 orbital, as a function of momentum transfer (q), and energy deposition (dE).
 Writes result to human-readable (and gnuplot-friendly) file, and/or binary.
 * For definitions/details, see: B.M.Roberts, V.A.Dzuba, V.V.Flambaum, M.Pospelov, Y.V.Stadnik,
 [Phys.Rev.D 93, 115037 (2016)](https://link.aps.org/doi/10.1103/PhysRevD.93.115037 "pay-walled");
 [arXiv:1604.04559](https://arxiv.org/abs/1604.04559 "free download").
 * Uses self-consistent Hartree Fock method
(optionally, can use parametric potential, which is faster but less accurate)
 * Note: need quite a dense grid [large number of points] for
   * a) highly oscillating J_L function at low r, and
   * b) to solve equation for high-energy continuum states.
 * Sums over 'all' continuum angular momentum states (and multipolarities)
   * Maximum values for l are input parameters

## dmeXSection

 * Calculates the cross-section and event rates for ionisation of atoms
 by scattering of DM particle.
 * Takes in the output of "atomicKernal"; is a "sub-program" of that work
 * Also calculates "observable" event rates, accounting for detector thresholds
 and resolutions. For now, just for DAMA detector. Will add for XENON

## wigner
 * Small routine to calculate 3,6,9-j symbols, and Clebsch Gordan coeficients
 * Either give input via command line directly (quote marks required)
   * e.g., _./wigner '<0.5 -0.5, 0.5 0.5| 1 0>'_
   * or e.g., _./wigner '(0.5 1 0.5, -0.5 0 0.5)'_ etc.
 * Or, give an input file, that contains any number of symbols, all on new line
   * e.g., _./wigner -f myInputFile.in_
   * nb: the '-f' flag can be dropped in the '.in' file extension is used
   * Do not use quote marks in input file. Lines marked '!' or '#' are comments
 * 3j symbols must start with '('; 6,9j with '{', and CG with '<' (this is how code knows which symbol to calculate).
 * but, each number can be seperated by any symbol (space, comma etc.)

## fitParametric

 * Finds the best-fit parameters for two-parameter parametric potentials
   (Green, or Tietz potentials)
 * Takes input/target states from fitParametric.in
