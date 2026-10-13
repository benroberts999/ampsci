\page module_matrixelements Matrix Elements

\brief Matrix elements of any operator, with optional RPA (core polarisation)

The `matrixElements` module (run as `Module::matrixElements{}`) calculates
reduced matrix elements (or hyperfine constants, or stretched-state matrix
elements) of any available operator, between valence
states, with optional RPA corrections.

* Runs for all pairs of valence states allowed by the selection rules
  (diagonal and/or off-diagonal).
* Optionally includes core states, and/or uses the spectrum (which includes
  correlations, if it was made with them) in place of the valence states.
* RPA may be included with the TDHF, diagram, or basis methods.
* Structure radiation and normalisation are _not_ included here: see
  \ref module_structurerad.

Frequency handling: frequency-dependent operators (e.g. E1v) are evaluated
at each pair's transition frequency, which is their physical frequency; the
`omega` option refers to the frequency RPA is solved at (`omega = each;`
re-solves RPA at each transition frequency). The rarely-needed
`omega_operator` option pins the operator at a fixed frequency instead.

Example:

```java
Module::matrixElements{
  operator = E1;
  rpa = TDHF;
  omega = 0.0;
}
```

The calculations are performed by the Amplitudes library: see the
documentation for \ref Amplitudes::matrix_elements "Amplitudes::matrix_elements()" and \ref Amplitudes::MEdata "Amplitudes::MEdata"
(the module, \ref Module::matrixElements "Module::matrixElements()", only parses input and prints).

Available options (from `./ampsci -m matrixElements`):

```java
// Available Module::matrixElements options/blocks
Module::matrixElements{
  // Matrix elements of any operator for HF/Brueckner
  // valence states. Supports RPA, diagonal and off-diagonal
  // elements, core states, and optional use of the spectrum
  // instead of valence states.
  // Note: SR and Norm are not included here; calculated
  // seperately in own modules
  operator;
    // e.g., E1, hfs (see ampsci -o for available operators)
  options{}
    // options specific to operator (see ampsci -o
    // 'operator')
  rpa;
    // Method used for RPA: true(=TDHF), false, TDHF, basis,
    // diagram [true]
  rpa_options{}
    // Block: some further options for RPA
  omega;
    // Text or number. Frequency RPA is solved at. Put
    // 'each' to solve at correct frequency for each
    // transition. [0.0]
  omega_operator;
    // Frequency-dependent operators are evaluated at the
    // transition frequency for each element, which is the
    // physical frequency. Set this to pin the operator at a
    // fixed frequency instead (rarely meaningful; mainly
    // for comparison to older calculations).
  printBoth;
    // print <a|h|b> and <b|h|a> [false]
  include_core;
    // If true, includes core states in calculation. Will
    // use HF core, unless use_spectrum is true [false]
  use_spectrum;
    // If true (and spectrum available), will use spectrum
    // for valence states [false]
  diagonal;
    // Calculate diagonal matrix elements (if non-zero)
    // [true]
  off-diagonal;
    // Calculate off-diagonal matrix elements (if non-zero)
    // [true]
  what;
    // What to calculate? Options are: Reduced (reduced
    // matric elements), Stetched (stretched states, with
    // j=m= [j=min(ja,jb) for off-diagonal]), or HFConstant
    // for (hyperfine A,B,etc. constants). Default is
    // Reduced, except for hyperfine operator, for which it
    // is HFConstant
}

// Available rpa_options options/blocks
rpa_options{
  eps;
    // Convergance goal [1.0e-10]
  eta;
    // Damping factor - be careful with this [0.35]
  max_iterations;
    // Maximum number of iterations. 1 should correspond to
    // first-order RPA [100]
}
```
