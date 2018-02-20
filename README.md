Dirac BS code.

Solves local central-field problem for Dirac equation.
Based on W. Johnson code, with a few extensions that improve numerical stability and accuracy.

With reasonable choices for the integration grids, typically converges to around 1e-16

The program named "h-like.cpp" is an example, that solves for H-like ions.
It takes a few input parameters from "h-like.in" file.
(For now, it is compiled with the shell script in h-like.b)

Includes a parameter to vary the effective speed of light - allowing non-relativistic approximation.

Note: makes use of GSL libraries: https://www.gnu.org/software/gsl/

[For example, with ubuntu: sudo apt-get install libgsl0-dev]
