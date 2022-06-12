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
  * Breit issue? (see below)

## Ladder diagrams (#13)
  * Form addition to Cor. Pot.
  * Screening?

## StrucRad + Diagram RPA
  * Option to use QkTable?
    * Can optionally calculate Qk table at Basis{}. Then, if we have it, use it
  * Store <a|h|b> and <a|dV|b> for each reqd pair?

## B-spline stability (#26)
  * Also: sign of some basis functions wrong sometimes?
  * Johnson-style splines not stable.. boundary conditions OK?
  * Better option for integration (form Hij)
  * Option to write Yk and/or Qk tables from splines block
    * Can then potentially run multiple basis blocks?

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
  * Use Qk,Pk functions? and/or QkTable [optionally]?
  * RPAD: minor "eps" race cond?

## Breit
  * More efficient (calc + store integrals, as for Yk?)
  * Breit matrix - for inclusion into Green's function (needs G?)

## Angular
  * Ck table a mess - merge with 6J?

## Double core Polarisation (#12)
  * Work with any class derivative?

## DiracODE (#24)
  * Write a general modern c++ DE solver
  * Allow non-local term (requires normalisation) (#11)
  * Efficient + numerically stable

## Coulomb
  * Make calling syntax consistent

## PNC
  * Work with diagram RPA
  * dV conj ???
   * Related to #3 ?
  * Different Sigma's

## Data format?
  * Standard data format for output/comparison?
  * Allow multiple runs
  * JSON?
  * Merge input files into input block?

## IO
  * Tidy + modernise
  * Write to disk: streamline
  * Binary: make class
  * Input block: to JSON?

## Unit tests
  * Improve + add more unit tests

## NIntegrate
  * Modernise + cleanup
