# To do

## Correlations (#21)

* Update class that stores, calculates sigma
* filename: label.sig2 -> CsI_label.sig2

## TDHF - physics

* Issue for even-parity operators (from de?) (#3)
  * Perhaps linked to (#11)
* dV conj - sometimes causes issues; + not consistent
* PNC: TDHF vs Diagram?
* Note: TDHF doesn't work for E1v, but diagram does!

## B-spline stability (#26)

* Also: sign of some basis functions wrong sometimes?
* Johnson-style splines not stable.. boundary conditions OK?
* Better option for integration (form Hij)
* Option to write Yk and/or Qk tables from splines block
  * Can then potentially run multiple basis blocks?

## Spectrum stability with Feynman Sigma (#28)

* Spectrum energies do not perfectly match with valence energies when using Feynman Sigma.
  * Seems related to hole-particle inclusion, but most likely numerical
* When using Goldstone sigma, however, they do.
* This implies it's probably a numerical error stemming from Feynman Sigma - but is not obvious.
* Typically difference is small, and doesn't seem to negatively affect much.
* However, when scaling sigma to exactly reproduce energy intervals, the scaling  is done for valence states, so the spectrum states will not match exactly!
* Possibly related to #22
* Possibly related to #26

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
* Faster Qk lookup
* Extra diagrams (triangle, 3rd-order)

## Diagram RPA

* Option to use QkTable? More consistent
  * Can optionally calculate Qk table at Basis{}. Then, if we have it, use it
* Store <a|h|b> and <a|dV|b> for each reqd pair?

## TDHF - performance

* Parallelised inefficiently (#6)
  * Due to thread-safe shared_ptr?
* Diagram RPA: use Qk table? Much slower?
* Write out to disk?
* Fewer allocations?
* Use Qk,Pk functions? and/or QkTable [optionally]?
* RPAD: minor "eps" race cond?

## Breit

* More efficient (calc + store integrals, as for Yk?)
* Breit matrix - for inclusion into Green's function (needs G?)
* Frequency-dependent Breit
* "Two-body" Breit (i.e., Breit into Sigma_2)

## Angular

* Ck table a mess - merge with 6J?

## Double core Polarisation (#12)

* Work with any class derivative?

## DiracODE (#24)

* Write a general modern c++ DE solver
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

## NIntegrate

* Modernise + cleanup
