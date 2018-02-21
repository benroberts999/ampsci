# Dirac bound state code.

Solves local central-field problem for Dirac equation.

Based on book by W. Johnson, with a few extensions that improve numerical
stability and accuracy.
 * W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007).


With reasonable choices for the integration grids, typically converges
to better than a few parts in 10^16

 * Includes an option to vary the effective speed of light -
allowing non-relativistic approximation.

Note: makes use of GSL libraries: https://www.gnu.org/software/gsl/

 * For example, with ubuntu: _sudo apt-get install libgsl0-dev_


## h-like.cpp

 * An example that solves for H-like ions
 * compile: _make h-like.x_
 * Input parameters/options in file: h-like.in


## fitParametric.cpp

 * Finds the best-fit parameters for two-parameter parametric potentials
 (Green, or Tietz potentials, See Johnson book]
 * Takes input/target states from fitParametric.in
 * compile: _make fitParametric.x_
