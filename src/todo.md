# To do

`[clang-format -i src/*.cpp src/*/*.cpp src/*/*.hpp]`

## Correlation
 * Construct correlation potential at specified energies
 * Easy way to use correct correlation potential matrix..?
    * e.g., in SolveMixedStates.....
 * Major cleanup..
 * filename: label.sig2 -> CsI_label.sig2
 * Option to use Sigma2 from Gold + rest from feyn

## Feynman
 * Polarisation operator instability: (#8)
   * Fails, e.g., for Cs with n_core < 3
 * Hole-particle; k=0 vs k=1?
 * Overall numerical stability
 * Exchange (2nd order): w1 and w1w2
 * g-part for Feynman
 * Breit issue?

## Ladder diagrams (#13)
 * Qk table
 * Full SDs coupled-cluster??
 * First-order ladder diagrams
   * What symmetries does l have?
 * Iterate ladder

## B-splines (#9)
 * Fix raw B-spline class
 * Fix B-spline basis:
 * Correct the r0 issue; correct number of splines
 * Implement Johnson version
 * Use more efficient integration?
   * Store coefs instead of expand?

## Hartree Fock
 * Re-write class (make less inter-dependent)
 * Option to use non-local DiracODE

## Wavefunction
 * Major cleanup
 * WF vs HF;

## TDHF
 * Issue for even-parity operators (from de?) (#3)
 * Perhaps linked to (#11)
 * Parallelised inefficiently (#6)
 * Diagram RPA: use Qk table?
 * Write out to disk?
 * Fewer allocations?
 * Use Qk,Pk functions?
 * RPAD: minor "eps" race cond?
 * PNC: TDHF vs Diagram?
 * dV conj - sometimes causes issues; + not consistent

## Double core Polarisation (#12)
 * Work with any class derivative?

## Operators (#20)
 * Fix up 'generate operator'
   * Have function: takes <userInputBlock> (+oper name)

## ADAMS
 * Write a general modern c++ DE solver
 * Allow non-local term (requires normalisation) (#11)
   * How to do for inward?
 * Efficient + numerically stable

## PNC
 * Work with diagram RPA
 * dV conj ???
 * Different Sigma's

## Continuum
 * Fix continuum class
 * Hartree-Fock?
 * Use splines?

## Data format?
  * Standard data format for output/comparison?
  * JSON?

## Overall
 * cleanup input options for main
 * Improve + add unit tests
 * CheckBlock: check for Blocks (not just options)
