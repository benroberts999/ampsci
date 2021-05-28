# To do

## Modules
 * Add g-factor (g-2) module

## Feynman
 * Polarisation operator instability: (#8)
   * Fails, e.g., for Cs with n_core < 3
 * Hole-particle; k=0 vs k=1?
 * Overall numerical stability
 * Exchange (2nd order): w1 and w1w2
 * Breit issue?


## Ladder diagrams (#13)
 * Qk table
 * Full SDs coupled-cluster??
 * First-order ladder diagrams
   * What symmetries does l have?
 * Iterate ladder

## B-splines (#9)
 * Correct the r0 issue; correct number of splines
 * Implement Johnson version
 * Use more efficient integration?
   * Store coefs instead of expand?

## TDHF
 * Issue for even-parity operators (from de?) (#3)
 * Parallelised inefficiently (#6)
 * Diagram RPA: use Qk table?
 * Write out to disk?
 * Fewer allocations?
 * Use Qk,Pk functions?
 * RPAD: minor "eps" race cond?

## Double core Polarisation (#12)
 * Work with any class derivative?


## Operators (#20)
 * Fix up 'generate operator'
   * Have function: takes <userInputBlock> (+oper name)

## Radiative potential (#10)
 * Re-write class in general form; make public
 * Allow different fittings for different l (or kappa)
 * Fix pointlike nucleus version (currently very slow)

## ADAMS
 * Write a general modern c++ DE solver
 * Efficient + numerically stable

## Continuum
 * Fix continuum class
 * Hartree-Fock?
 * Use splines?

## Data format?
  * Standard data format for output/comparison?
  * JSON?

## Overall
 * cleanup input options for main
 * Clean Wavefunction class
 * Clean DiracSpinor class
 * Improve + add unit tests
 * CheckBlock: check for Blocks (not just options)
