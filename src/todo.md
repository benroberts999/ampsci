# To do

## Correlations (#21)
  * Construct correlation potential at specified energies
  * Easy way to use correct correlation potential matrix..?
    * e.g., in SolveMixedStates.....
  * Major cleanup..
  * filename: label.sig2 -> CsI_label.sig2
  * Option to use Sigma2 from Gold + rest from feyn
    * Method = {Feynman, Goldstone, Mixed}
    * Feynman: include scr+hp by default
      * Calc fk if not given?

## Feynman (#22)
  * Cleanup, and re-do tests
  * Polarisation operator instability: (#8)
    * Fails, e.g., for Cs with n_core < 3
  * Hole-particle; k=0 vs k=1?
  * Overall numerical stability
  * Exchange (2nd order): w1 and w1w2
  * g-part for Feynman?
  * Breit issue?

## Ladder diagrams (#13)
  * Check eqs; try rho.
  * Form addition to Cor. Pot.

## StrucRad + Diagram RPA
  * Option to use QkTable?
  * Store <a|h|b> and <a|dV|b> for each reqd pair?

## B-spline stability (#26)
  * Also: sign of some basis functions wrong sometimes?
  * Johnson-style splines not stable.. boundary conditions OK?
  * Better option for integration (form Hij)
  * Option to write Yk and/or Qk tables from splines block
    * Can then potentially run multiple basis blocks?

## Breit
  * Very slow; use Yk table?
  * Form matrix (for Feyn) = needs G?

## Hartree Fock (#23)
  * Re-write class (make less inter-dependent)
  * Option to use non-local DiracODE (#11)
  * Two versions of Vdir (one in wf, one in HF) - BAD.

## Wavefunction (#23)
  * Major cleanup
  * WF vs HF;

## TDHF - physics
  * Mixed-states work with non-local DiracODE? (#11)
  * Issue for even-parity operators (from de?) (#3)
  * Perhaps linked to (#11)
  * dV conj - sometimes causes issues; + not consistent
  * PNC: TDHF vs Diagram?
  * Note: TDHF doesn't work for E1v, but diagram does!

## TDHF - performance
  * Parallelised inefficiently (#6)
  * Diagram RPA: use Qk table?
  * Write out to disk?
  * Fewer allocations?
  * Use Qk,Pk functions?
  * RPAD: minor "eps" race cond?

## Double core Polarisation (#12)
  * Work with any class derivative?

## DiracOperator (#20)
  * Fix up 'generate operator'
    * Have function: takes <userInputBlock> (+oper name)

## DiracODE (#24)
  * Write a general modern c++ DE solver
  * Allow non-local term (requires normalisation) (#11)
  * Efficient + numerically stable

## PNC
  * Work with diagram RPA
  * dV conj ???
  * Different Sigma's

## Data format?
  * Standard data format for output/comparison?
  * JSON?
  * Merge input files into input block?

## IO
  * Tidy
  * Write to disk: streamline
  * Binary: make class

## Overall
  * cleanup input options for main
  * Improve + add unit tests

## ampsci (ampsci.cpp)
  * Clean/tidy: separate input blocks
