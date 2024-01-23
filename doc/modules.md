# Modules

:: Descpription of modules system: available modules and options

- The modules system allows the easy calculation of any atomic properties after the wavefunction has been calculated.
- Any number of _modules_ can be run by adding `ModulebriefmoduleName{}' blocks.
- The code is designed so that you can easily create your own modules. See [doc/writing_modules.md](/doc/writing_modules.md) for details

Get a list of available modules: `./ampsci -m`

They are (at time of writing):

Available modules:

- Tests
  - Some basic wavefunction numerical tests
- testBasis
  - Tests of basis and spectrum
- WriteOrbitals
  - Write orbitals to disk for plotting
- matrixElements
  - Calculates matrix elements of any operator
- CI_matrixElements
  - Calculates matrix elements of any operator for CI wavefunctions
- lifetimes
  - Calculate radiative lifetimes (E1, E2, M1)
- polarisability
  - Calculates static polarisabilities
- dynamicPolarisability
  - Calculates dynamic polarisabilities
- transitionPolarisability
  - Calculates transition polarisabilities
- structureRad
  - Calculates Struct. Rad + Normalisation corrections to MEs
- fieldShift
  - Calculates field-shift constants (isotope shift)
- QED
  - QED corrections to energies/matrix elements
- Breit
  - Breit corrections to energies
- ladder
  - Calculates ladder diagrams and energy corrections
- Kionisation
  - Calculate atomic ionisation form-factors
- continuum
  - Compute and use continuum wavefunctions
- HFAnomaly
  - Calculates Bohr-Weisskopf effect and hyperfine anomaly
- screeningFactors
  - Calculates Feynman electron screening factors
- pnc
  - Calculates APV amplitudes
- exampleModule
  - A short description of the module

You can see all the available options by setting the 'help' option, e.g.,

- `ModulebriefmoduleName{ help; }`

You can also get most of this information directly from the command-line:

- `./ampsci -m  <ModuleName>`
  - Prints list of available Modules (same as --modules)
  - ModuleName is optional. If given, will list available options for that Module
