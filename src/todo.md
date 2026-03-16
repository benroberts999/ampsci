# To do

## CI

* Allow use of Yk, instead of Qk, to save memory
* Update solutions:
  * re-calculate more?
* Option for E-denominator type
* Screening for Sigma(2)
  * fk screening
  * etak hp

## TDHF

* Issue for even-parity operators (from de?) (#3)
  * Green's method fails for even operators when k=k'
* Also: check correctness in general case (#33)
  * Note: TDHF doesn't work for E1v, but diagram does!

## Qk, Coulomb, Angular

* QkTable - better hash table?
* Currently, lookup is a major bottle-neck

## Feynman (+correlations)

* Cleanup, and re-do tests
* Polarisation operator instability: (#8)
  * Fails, e.g., for Cs with n_core < 3
* Overall numerical stability
* Exchange (2nd order): w1 and w1w2 (#3)
* Breit issue? (#34)

## Breit

* Frequency-dependent Breit
* "Two-body" Breit into Sigma(2) - implemented, but needs testing (#31)

## Modules

* Important modules (e.g., PNC): clean + improve
* System for incorperating external modules?
* Compile seperately: required read wf from disk!
  * and/or, call external functions

## Ladder diagrams (#36)

* Implement and test ladder corrections
* Extra diagrams (triangle, 3rd-order)
  * Complete to third-order?

## Double-core polarisation

* To PNC, but also to polarisability

## Operators

* Add important operators (NSD-PNC, Yakawa, LVI)
* Possibility for generalised 1 and 2 particle operators?
* Core pol should work like operator!

## Performance

* TDHF
  * Parallelised inefficiently (#6): thread-safe shared_ptr?
  * Fewer allocations?
  * Particularly noticable with Breit
* Struc Rad
  * Inneficient, not parellised efficiently?

## Improve/modernise

* Consider TOML input? JSON?
* Clean UserInput
* Standardise input/output
* Option for JSON output?
