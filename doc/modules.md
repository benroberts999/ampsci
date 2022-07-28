# Modules

\brief Descpription of modules system: available modules and options

[[Home](/README.md)]

- The modules system allows the easy calculation of any atomic properties after the wavefunction has been calculated.
- Any number of _modules_ can be run by adding `Module::moduleName{}' blocks.
- The code is designed so that you can easily create your own modules. See [doc/writing_modules.md](/doc/writing_modules.md) for details

Get a list of available modules: `$ ./ampsci -m`

They are:

- Tests
- WriteOrbitals
- AtomicKernal
- BohrWeisskopf
- HFAnomaly
- HF_rmag
- screeningFactors
- BW_eta_sp
- pnc
- vertexQED
- QED
- testFeynman
- matrixElements
- lifetimes
- polarisability
- dynamicPolarisability
- transitionPolarisability
- structureRad
- fieldShift
- continuum
- ladder

You can see all the available options by setting the 'help' option, e.g.,

- `Module::moduleName{ help; }`

You can also get most of this information directly from the command-line:

- `$ ./ampsci -m  <ModuleName>`
  - Prints list of available Modules (same as --modules)
  - ModuleName is optional. If given, will list avaiable options for that Module

Details for each module are given below.
It's usually better to use the 'help' option of the code, as it will always be up-to-date, while this document may be out of date.

--------------------------------------------------------------------------------

## Module::matrixElements

Calculates matrix elements of the given operator

```cpp
Module::matrixElements{
  operator; // e.g., E1, hfs
  options; // options specific to operator; blank by dflt
  rpa; // true(=TDHF), false, TDHF, basis, diagram
  omega; // freq. for RPA
  radialIntegral; // false by dflt (means red. ME)
  printBoth; // print <a|h|b> and <b|h|a> (dflt false)
  onlyDiagonal; // only <a|h|a> (dflt false)
}
```

- `$./ampsci -o` will print a list of available operators
- Some operators take further options. Set `options{help;}` to see list
- `$ ./ampsci -o <OperatorName>`
  - Prints list of available operators (same as --operators)
  - OperatorName is optional. If given, will list avaiable options for Operator

--------------------------------------------------------------------------------

## Module::lifetimes

Lifetimes:
Note: Uses _valence_ states - so, must ensure all lower states have been included in the valence list for accurate results.

Available Module::lifetimes

```cpp
Module::lifetimes{
  E1; // Include E1 transitions? [true]
  E2; // Include E2 transitions? [false]
  rpa; // Include RPA? [true]
  StrucRadNorm; // Include SR+Norm correction (only for E1)? [false]
}
```

--------------------------------------------------------------------------------

## Module::polarisability

Calculate atomic polarisabilities at single frequency

```cpp
Module::polarisability{
  rpa; // Include RPA? [true]
  omega; // frequency (for single w) [0.0]
  tensor; // Also calculate tensor alpha_2(w) (as well as a0) [false]
  drop_continuum; // Discard states from the spectrum with e>0 - these can cause spurious resonances [false]
  drop_states; // List. Discard these states from the spectrum for sum-over-states for valence part of alpha, and from TDHF by orthogonality (must be in core/valence) []
  StrucRad; // SR: include SR+Norm correction [false]
  n_min_core; // SR: Minimum n to include in SR+N [1]
  max_n_SR; // SR: Maximum n to include in the sum-over-states for SR+N [9]
  Qk_file; // SR: filename for QkTable file. If blank will not use QkTable; if exists, will read it in; if doesn't exist, will create it and write to disk. Save time (10x) at cost of memory.
}
```

--------------------------------------------------------------------------------

## Module::dynamicPolarisability

Calculate atomic dynamic polarisabilities

```cpp
Module::dynamicPolarisability{
  tensor; // Do tensor polarisability a2(w) (as well as a0) [false]
  rpa; // Include RPA? [true]
  core_omega; // Frequency-dependent core? If true, core part evaluated at each frequency. If false, core evaluated once at w=0 [true]
  rpa_omega; // Frequency-dependent RPA? If true, RPA solved at each frequency. If false, RPA solved once at w=0 [true]
  num_steps; // number of steps for dynamic polarisability [10]
  omega_minmax; // list (frequencies): omega_min, omega_max (in au) [0.01, 0.1]
  lambda_minmax; // list (wavelengths, will override omega_minmax): lambda_min, lambda_max (in nm) [600, 1800]
  method; // Method used for dynamic pol. for a0(w). Either 'SOS' (sum-over-states) or 'MS' (mixed-states=TDHF). MS can be unstable for dynamic pol. [SOS]
  replace_w_valence; // Replace corresponding spectrum states with valence states - circumvents spectrum issue! [false]
  drop_continuum; // Discard states from the spectrum with e>0 - these can cause spurious resonances [false]
  drop_states; // List. Discard these states from the spectrum for sum-over-states []
  filename; // output filename for dynamic polarisability (a0_ and/or a2_will be appended to start of filename) [identity.txt (e.g., CsI.txt)]
  StrucRad; // SR: include SR+Norm correction [false]
  n_min_core; // SR: Minimum n to include in SR+N [1]
  max_n_SR; // SR: Maximum n to include in the sum-over-states for SR+N [9]
  Qk_file; // SR: filename for QkTable file. If blank will not use QkTable; if exists, will read it in; if doesn't exist, will create it and write to disk. Save time (10x) at cost of memory.
}
```

--------------------------------------------------------------------------------

## Module::transitionPolarisability

Calculate transition polarisabilities (just alpha at the moment)

```cpp
Module::transitionPolarisability{
  transition; // List. states (e.g., 6s,6s) []
  rpa; // Include RPA? [true]
  omega; // frequency (for single w) [default: transition freq.]
  StrucRad; // SR: include SR+Norm correction [false]
  n_min_core; // SR: Minimum n to include in SR+N [1]
  max_n_SR; // SR: Maximum n to include in the sum-over-states for SR+N [9]
  Qk_file; // SR: filename for QkTable file. If blank will not use QkTable; if exists, will read it in; if doesn't exist, will create it and write to disk. Save time (10x) at cost of memory.
}
```

--------------------------------------------------------------------------------

## Module::structureRad

Calculates Structure Radiation + Normalisation corrections to matrix elements

```cpp
Module::structureRad{
  operator; // e.g., E1, hfs
  options; // options specific to operator; blank by dflt
  rpa; // true(=TDHF), false, TDHF, basis, diagram
  omega; // freq. for RPA
  printBoth; // print <a|h|b> and <b|h|a> (dflt false)
  onlyDiagonal; // only <a|h|a> (dflt false)
  Qk_file; // filename for QkTable file. If blank will not use QkTable; if exists, will read it in; if doesn't exist, will create it and write to disk. Save time (10x) at cost of memory. Note: Using QkTable implies splineLegs=true
  n_minmax; // list; min,max n for core/excited: (1,inf)dflt
  splineLegs; // Use splines for diagram legs (false dflt)
}
```

--------------------------------------------------------------------------------

## Module::fieldShift

Calculates field-shift factor, F = d(E)/d(<r^2>)

```cpp
Module::fieldShift{
  /*
  Calculates field shift: F = d(E)/d(<r^2>)
  */
  print; // Print each step? [true]
  min_pc; // Minimum percentage shift in r [1.0e-3]
  max_pc; // Maximum percentage shift in r [1.0]
  num_steps; // Number of steps for derivative (for each sign)? [10]
}
```

--------------------------------------------------------------------------------

## Module::Tests

Module to perform sum run-time tests

```cpp
Module::Tests{
  orthonormal; // 
  orthonormal_all; // 
  Hamiltonian; // 
  boundaries; // 
  basis; // 
}
```

--------------------------------------------------------------------------------

## Module::WriteOrbitals

Writes orbitals to a text file for plotting (gnuplot format)

```cpp
Module::WriteOrbitals{
  label; //
  l; //
}
```

--------------------------------------------------------------------------------

## Module::BohrWeisskopf

Calculates BW effect
WARNING: This may need to be re-checked after refactor!

```cpp
Module::BohrWeisskopf{
  rpa; //
  rpa_diagram; //
  screening; //
  hfs_options; //
}
```

--------------------------------------------------------------------------------

## Module::HFAnomaly

Calculates HF Anomaly effect

```cpp
Module::HFAnomaly{
  rpa; //
  options; //
  A; //
}
```

--------------------------------------------------------------------------------

## Module::HF_rmag

Tuning Rmag to fit hyperfine anomaly

```cpp
Module::HF_rmag{
  n; //
  kappa; //
  A2; //
  1D2; //
  rpa; //
  num_steps; //
  mu1; //
  mu2; //
  I1; //
  I2; //
  eps_targ; //
  e1; //
  e2; //
}
```

--------------------------------------------------------------------------------

## Module::BW_eta_sp

Calculates ets_sp function for transating s-p BW effect

```cpp
Module::BW_eta_sp{
}
```

--------------------------------------------------------------------------------

## Module::pnc

Calculates PNC amplitudes

```cpp
Module::pnc{
  t; //
  c; //
  transition; //
  nmain; //
  rpa; //
  omega; //
  E1_rpa_it; //
  pnc_rpa_it; //
}
```

--------------------------------------------------------------------------------

## Module::vertexQED

In progress: calculates QED vertex correction using approx formula

```cpp
Module::vertexQED{
  operator; // operator (e.g., E1 or hfs)
  options; // operator options (same as matrixElements)
  rrms; // nuclear rms, for QED part
  onlyDiagonal; // only print <a|h|a>
  radialIntegral; // false by default (means red. mat. el)
  A_vertex; // A vtx factor; blank=default
  b_vertex; // A vtx factor; =1 by default
  rpa; // include RPA? NOT USED FOR NOW
  omega; // freq. for RPA; NOT USED FOR NOW
}
```

--------------------------------------------------------------------------------

## Module::QED

Calculates QED corrections to energies, matrix elements

```cpp
Module::QED{
  A_vertex; // A vtx factor; blank means dflt
  b_vertex; // B vtx factor; =1 by default
  rrms; // double; effective rrms used in radiative potential
  out_file; // Results appended to this file (if given)
  matrixElements; // sub-block; same as Module::matrixElements
  coreQED; // bool; Include QED into core? Or just valence
  scale_l; // list; Scale factors for each l; e.g., 1,0 means s, but no p,d. default = 1, meaning include all
  vertex; // bool; Calculate vertex corrections?
}
```

--------------------------------------------------------------------------------

## Module::continuum

Calculates continuum wavefunctions, optionally matrix elements, and writes wfs to disk

```cpp
Module::continuum{
  energy; // List. energy for cntm states (>0) [0.5]
  max_l; // maximum l
  filename; // filename for output
  operator; // Operator to calculate matrix elements (e.g., E1)
  options; // options specific to operator; blank by dflt
  rpa; // Include RPA (TDHF for now)? [false]
  omega; // Frequency for RPA [0.0]
}
```

--------------------------------------------------------------------------------

## Module::ladder

Calculates ladder diagrams, stores in file (Lfile).

```cpp
Module::ladder{
  min; // lowest core n to include [0]
  max; // maximum excited n to include [99]
  max_l; // maximum excited l to include [99]
  max_k; // maximum k to include in Qk [99]
  include_L4; // Inlcude 4th Ladder diagram [false]
  fk; // List of doubles. Effective screening factors. Used to calculate Lk. []
  eta; // List of doubles. Effective hp factors. Only used to print energy shift. []
  Qfile; // filename to read/write Qk integrals
  form_Q; // Form or read Qk? (if have lk already, dont' need!) [true]
  Lfile; // filename to read/write Qk integrals
  progbar; // Print progress bar? [true]
  max_it; // Max # iterations [15]
  eps_target; // Target for convergance [1.0e-4]
}
```
