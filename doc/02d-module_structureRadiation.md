\page module_structurerad Structure Radiation

\brief Structure radiation, normalisation of states, and Brueckner orbital corrections to matrix elements

The `StructureRadiation` module (run as `Module::StructureRadiation{}`)
calculates the second-order (in the residual Coulomb interaction) MBPT
corrections to matrix elements:

* Structure radiation (SR): the "top", "bottom", and "centre" diagrams
* Normalisation of states
* Brueckner orbital (BO) corrections (included only when the external legs
  are not already Brueckner orbitals)

The total corrected matrix element is
\f[
  t^{\rm tot}_{ab} = t^{(0)}_{ab} + \delta V_{ab} + \delta t^{\rm SR}_{ab}
    + \delta t^{\rm Norm}_{ab} + \delta t^{\rm BO}_{ab}.
\f]

Aliases: the module may also be run as `Module::StructureRad{}` or
`Module::StrucRad{}`.

* Requires a basis (B-splines).
* The `legs` option selects which states are used for the external lines
  (HF valence, basis, spectrum, or Brueckner orbitals).
* The `Qk_file` option caches the Coulomb integrals to disk: roughly 10x
  faster, at the cost of memory. An existing file is re-used, and any
  missing integrals are calculated and appended automatically.
* Frequency handling is as for the \ref module_matrixelements module.
  Frequency-dependent operators are at each pair's transition frequency
  (`omega_operator` pins them); `omega` refers to the RPA and the SR
  tables/denominators, whose frequency dependence is usually weak.

Example:

```java
Module::StructureRadiation{
  operator = hfs;
  rpa = false;
  Qk_file = true;
}
```

The calculations are performed by \ref Amplitudes::matrix_element "Amplitudes::matrix_element()" (see
also \ref Amplitudes::MEdata "Amplitudes::MEdata"), using
\ref MBPT::StructureRad "MBPT::StructureRad"; the module only parses input
and prints. See \ref MBPT::StructureRad "MBPT::StructureRad" for the
diagram formulas [Johnson et al., At. Data Nucl. Data Tables 64, 279
(1996)].

Available options (from `./ampsci -m StructureRadiation`):

```java
// Available Module::StructureRadiation options/blocks
Module::StructureRadiation{
  // Calculates structure radiation, normalisation of
  // states, and Brueckner orbital corrections to matrix
  // elements using perturbation theory
  operator;
    // e.g., E1, hfs
  options{}
    // options specific to operator; blank by dflt
  rpa;
    // true(=TDHF), false, TDHF, basis, diagram [true]
  omega;
    // Text or number. Frequency for RPA and the SR
    // tables/denominators. Put 'each' to solve at correct
    // frequency for each transition. [0.0]
  omega_operator;
    // Frequency-dependent operators are evaluated at the
    // transition frequency for each element, which is the
    // physical frequency. Set this to pin the operator at a
    // fixed frequency instead (rarely meaningful; mainly
    // for comparison to older calculations).
  printBoth;
    // print <a|h|b> and <b|h|a> (dflt false)
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
  Qk_file;
    // true/false/filename - SR: filename for QkTable file.
    // If blank will not use QkTable; if exists, will read
    // it in; if doesn't exist, will create it and write to
    // disk. If 'true' will use default filename. Save time
    // (10x) at cost of memory. Note: Using QkTable implies
    // legs=basis
  n_minmax;
    // list; min,max n for core/excited: (1,inf)dflt
  k_cut;
    // Maximum multipolarity k to include in Coulomb Qk.
    // Default: all
  include_core;
    // If true, includes core states in calculation. Will
    // use HF core, unless legs=spectrum [false]
  legs;
    // Which states to use for diagram legs?
    // hf/basis/spectrum/brueckner [hf]
  fk;
    // list: Coulomb screening factors for each k
  etak;
    // list: Hole-particle factors for each k
}
```
