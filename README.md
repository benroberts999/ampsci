# Relativistic, self-consistent atomic structure code.

Solves local central-field problem for Dirac equation, and does Hartree method.
(For now, does not include exchange).

With reasonable choices for the integration grids, typically converges
to better than a few parts in 10^16

The part that solves the Dirac eigenvalue DE is based on book by W. Johnson,
with a few extensions that improve numerical stability and accuracy.
 * W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007).

 * Includes an option to vary the effective speed of light -
allowing non-relativistic approximation.

 * Wavefunctions are in form psi = (1/r) [iP,Q], (using Dirac basis)

 * Also, has ability to solve for continuum states

Note: makes use of GSL libraries: https://www.gnu.org/software/gsl/

 * For example, with ubuntu: _$sudo apt-get install libgsl0-dev_
 * Probably also needs LAPACK/BLAS: _$sudo apt-get install libblas-dev libatlas-dev liblapack-dev_


### Compiling and use:

 * All programs compiled using the Makefile
 (run _$make_ or _$make programName.x_)
 * All executables end with '.x' suffix; run like _$./programName.x_
 * All programs have input options, stored and read from 'programName.in' file

## h-like.x

 * An example that solves for H-like ions

## hartree.x

 * Solves Hartree potential (no exchange) for core + valence states
 * Takes core configuration: Noble gas + extra. end with ' . ' (needs space)
 E.g.:
   * For Cs: 'Xe . '
   * For Au: 'Xe 4f14 5d10 . '
   * For Tl: 'Xe 4f14 5d10 6s2 . '
 * Solve single-electron valence states in the Hartree potential of given core
   * For now, the core must be closed-shell
 * Includes finite nuclear size (assumes spherical nucleus)
 * As of yet, does not write wavefunctions to disk
 * NOTE: Not 100% finished yet, still has some issues. Works well enough though

## atomicKernal.x

 * Calculates the "Atomic Kernal" (for scattering/ionisation) for each core
 orbital, as a function of momentum transfer (q), and energy deposition (dE).
 Writes result to human-readable (and gnuplot-friendly) file.
 * For definitions/details, see: B.M.Roberts, V.A.Dzuba, V.V.Flambaum, M.Pospelov, Y.V.Stadnik, 
 [Phys.Rev.D 93, 115037 (2016)](https://link.aps.org/doi/10.1103/PhysRevD.93.115037 "pay-walled"); 
 [arXiv:1604.04559](https://arxiv.org/abs/1604.04559 "free download").
 * Uses self-consistent Hartree method (optionally, can use parametric potential, which is faster but less accurate)
 * Note: need quite a dense grid [large number of points] for
   * a) highly oscillating J_L function at low r, and
   * b) to solve equation for high-energy continuum states.
   * NGP=25000 allows good dE and q range, and is fairly quick
 * Sums over 'all' continuum angular momentum states (and multipolarities)
   * Maximum values for l are input parameters
 * Yes, I know I spelled "kernel" wrong..

## parametricPotential.x

 * Solves Dirac equation using Green/Tietz parametric potentials
 * You can give it parameters (H,g,t,d), or it will use defaults
 * Optionally: give it the core configuration (As in 'hartree' program)

## fitParametric.x

 * Finds the best-fit parameters for two-parameter parametric potentials
   (Green, or Tietz potentials, See Johnson book)
 * Takes input/target states from fitParametric.in
