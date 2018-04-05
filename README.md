# Dirac bound state code.

Solves local central-field problem for Dirac equation.

Based on book by W. Johnson, with a few extensions that improve numerical
stability and accuracy.
 * W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007).


With reasonable choices for the integration grids, typically converges
to better than a few parts in 10^16

 * Includes an option to vary the effective speed of light -
allowing non-relativistic approximation.

 * Also, has ability to solve for continuum states

Note: makes use of GSL libraries: https://www.gnu.org/software/gsl/

 * For example, with ubuntu: _sudo apt-get install libgsl0-dev_

### Compiling and use:

 * All programs compiled using the Makefile (run $make or $make programName.x)
 * All exectables end with '.x' suffix
 * All programs have input options, stored and read from 'programName.in' file

## h-like.x

 * An example that solves for H-like ions
 * compile: _make h-like.x_
 * Input parameters/options in file: h-like.in

## fitParametric.x

 * Finds the best-fit parameters for two-parameter parametric potentials
 (Green, or Tietz potentials, See Johnson book]
 * Takes input/target states from fitParametric.in
 * compile: _make fitParametric.x_

## parametricPotential.x

 * Solves Dirac equation using Green/Tietz parametric potentials
 * You can give it parameters (H,g,t,d), or it will use defaults
 * Optionally: give it the core configuration Eg. for Ag, core is:
   Kr 4d10 .
