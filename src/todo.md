# To do

## Correlations (#21)

* Update class that stores, calculates Sigma
* Option to mix methods (e.g., Goldstone for g-part?)

## TDHF - physics (#3)

* Issue for even-parity operators (from de?) (#3)
  * Green's method fails for even operators when k=k'
* Also: check correctness in general case
  * Note: TDHF doesn't work for E1v, but diagram does!

## TDHF - code

* Streamline the core-polarisation classes
* Take unified bunch of options

## Qk and Coulomb

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
* Hole-particle; k=0 vs k=1?
* Overall numerical stability
* Exchange (2nd order): w1 and w1w2
* g-part for Feynman?
* Breit issue? (see below)

## Breit

* Frequency-dependent Breit
* "Two-body" Breit into Sigma(2)
* Breit matrix - for inclusion into Green's function (needs G?)

## Modules

* PNC
* Cleanup and improve other modules
* System for incorperating external modules?
* HF_anomaly - incorperate into single + clean

## Ladder diagrams (#13)

* Implement and test ladder corrections
* Extra diagrams (triangle, 3rd-order)
  * Complete to third-order?

## Double core Polarisation (#12)

* Work with any class derivative?

## Performance

* TDHF
  * Parallelised inefficiently (#6): thread-safe shared_ptr?
  * Fewer allocations?
  * Particularly noticable with Breit
* Breit: ineficient?
* Qk table: better table?

## Improve/modernise

* Standardise input/output
* Option for JSON output?
* NIntegrate
