# Modules

:::: Descpription of modules system: available modules and options

- The modules system allows the easy calculation of any atomic properties after the wavefunction has been calculated.
- Any number of _modules_ can be run by adding `Module::moduleName{}' blocks.
- The code is designed so that you can easily create your own modules. See [doc/writing_modules.md](/doc/writing_modules.md) for details

Get a list of available modules: `./ampsci -m`

Output will look something like this:

```txt
Available modules: 
 * Tests
     Some basic wavefunction numerical tests
 * testBasis
     Tests of basis and spectrum
 * WriteOrbitals
     Write orbitals to disk for plotting
 * matrixElements
     Calculates matrix elements of any operator
 * CI_matrixElements
     Calculates matrix elements of any operator for CI wavefunctions
 * thirdOrderME
     Calculates Third-order matrix elements
 * lifetimes
     Calculate radiative lifetimes (E1, E2, M1)
 * polarisability
     Calculates static polarisabilities
 * dynamicPolarisability
     Calculates dynamic polarisabilities
 * transitionPolarisability
     Calculates transition polarisabilities
 * structureRad
     Calculates Struct. Rad + Normalisation corrections to MEs
 * fieldShift
     Calculates field-shift constants (isotope shift)
 * QED
     QED corrections to energies/matrix elements
 * Breit
     Breit corrections to energies
 * ladder
     Calculates ladder diagrams and energy corrections
 * Kionisation
     Calculate atomic ionisation form-factors
 * continuum
     Compute and use continuum wavefunctions
 * HFAnomaly
     Calculates Bohr-Weisskopf effect and hyperfine anomaly
 * screeningFactors
     Calculates Feynman electron screening factors
 * pnc
     Calculates APV amplitudes
 * muonPV
     For testing/playing with muonic PV
 * VQE
     For testing/playing with VQE method
 * exampleModule
     A short description of the module
```

You can also get most of this information directly from the command-line:

- `./ampsci -m  <ModuleName>`
  - Prints list of available Modules (same as --modules)
  - ModuleName is optional. If given, will list available options for that Module
  - Note the output is in the same format as required by the input file - you can copy+paste this into your input file.

```sh
./ampsci -m MatrixElements
```

```java
// Available Module::MatrixElements options/blocks
Module::MatrixElements{
  // e.g., E1, hfs (see ampsci -o for available operators)
  operator;
  // options specific to operator (see ampsci -o 'operator')
  options{}
  // Method used for RPA: true(=TDHF), false, TDHF, basis, diagram [true]
  rpa;
  // Text or number. Freq. for RPA (and freq. dependent operators). Put 'each'
  // to solve at correct frequency for each transition. [0.0]
  omega;
  // print <a|h|b> and <b|h|a> [false]
  printBoth;
  // If true (and spectrum available), will use spectrum for valence states
  // [false]
  use_spectrum;
  // Calculate diagonal matrix elements (if non-zero) [true]
  diagonal;
  // Calculate off-diagonal matrix elements (if non-zero) [true]
  off-diagonal;
  // Options for Structure Radiation and normalisation (details below)
  StructureRadiation{}
}


// Available StructureRadiation options/blocks
StructureRadiation{
  // If this block is included, SR + Normalisation corrections will be included

  // true/false/filename - SR: filename for QkTable file. If blank will not use
  // QkTable; if exists, will read it in; if doesn't exist, will create it and
  // write to disk. If 'true' will use default filename. Save time (10x) at cost
  // of memory. Note: Using QkTable implies splines used for diagram legs
  Qk_file;
  // list; min,max n for core/excited: [1,inf]
  n_minmax;
}
```
