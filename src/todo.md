# To do

## CI

* Update solutions:
  * re-calculate more?
  * Check if any didn't converge
  * Control over tolerance?
* Option for E-denominator type
* Screening for Sigma(2)
  * fk screening
  * etak hp
* Breit into g_2

## Correlations (#21)

* Update full integration tests
* Clean up interface (ampsci/wavefunction)

## TDHF - physics (#3)

* Issue for even-parity operators (from de?) (#3)
  * Green's method fails for even operators when k=k'
* Also: check correctness in general case
  * Note: TDHF doesn't work for E1v, but diagram does!

## CorePolarisation - code

* Streamline the core-polarisation classes: unified calling syntax
* Allow it to work like an operator!

## Qk, Coulomb, Angular

* QkTable - better hash table?
* Create QkTable at Basis{}?
  * Option to use different basis for RPA/Sigma etc?

## B-spline stability (#26)

* Improve spline stability:
  * Sign of some basis functions wrong sometimes?
  * Johnson-style splines not stable.. boundary conditions OK?
  * Better option for integration (form Hij)?
* Spectrum stability (#28)
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
* Overall numerical stability
* Exchange (2nd order): w1 and w1w2
* g-part for Feynman?
* Breit issue? (see below)
* New method for Green's function at imag. energy
* Mixed methods: etak and fk

## Breit

* Frequency-dependent Breit
* "Two-body" Breit into Sigma(2)
* Breit matrix - for inclusion into Green's function (needs G?)
* Breit into RPA(D)/RPA - test?

## Modules

* PNC
* Cleanup and improve other modules
* System for incorperating external modules?
* Compile seperately: required read wf from disk!
  * and/or, call external functions
* Compile seperate modules: use ampsci as a library

## Ladder diagrams (#13)

* Implement and test ladder corrections
* Extra diagrams (triangle, 3rd-order)
  * Complete to third-order?

## Performance

* TDHF
  * Parallelised inefficiently (#6): thread-safe shared_ptr?
  * Fewer allocations?
  * Particularly noticable with Breit

## Improve/modernise

* Consider TOML input?
* Clean UserInput
* Standardise input/output
* Option for JSON output?
