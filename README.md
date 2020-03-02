# Relativistic, self-consistent atomic structure code.

Solves the Dirac equation for atomic systems using the Hartree-Fock method.
Fully relativistic, includes finite-nuclear size, and can solve for continuum states (energy normalisation).


### Compilation:

 * All programs compiled using the Makefile (run _$make_ or _$make programName_)
 * Install make on ubutnu: _$sudo apt-get install make_
 * Tested with g++ and clang++ on linux and mac (requires c++17)

Note: makes use of GSL libraries (tested with ver:2.4): https://www.gnu.org/software/gsl/, and LAPACK. These must be installed for the code to run.

  * For example, install GSL with ubuntu: _$sudo apt-get install libgsl-dev_
  * And LAPACK/BLAS: _$sudo apt-get install libatlas-base-dev liblapack-dev libblas-dev_
  * On mac: Just use homebrew to install gsl: _$brew install gsl_
  * (homebrew is a package manager; install from https://brew.sh/)

For windows, the easiest way (for me, anyway) is to make use of the recent 'windows subsystem for linux'. Instructions on installation/use here: https://www.roberts999.com/posts/2018/11/wsl-coding-windows-ubuntu
Then, the compilation + use can proceed as per above.

 * **NOTE:** If you get the following error message on compile:
_error: unsupported option -fopenmp_
change '_UseOpenMP=yes_' to '_UseOpenMP=no_' in Makefile

 * Sometimes, the compiler will not be able to find the correct libraries (particular, e.g., on clusters). In this case, there are two options in the Makfefile: **ExtraInclude** and **ExtraLink**
 * These add paths to include the correct directories for both -I "includes" (for compilation), and -L link flags (for linking libraries) in Makefile. These can be a little tricky to get right (don't include the -I or -L)
 * The current defaults are setup to link/compile correctly on UQ's getafix server - you just need to uncomment the two lines in Makefile (remember to load the correct getafix modules, we need 'gnu' and 'gsl')
 * If you get this error (for example), that's why:
  _error: too few arguments to function ‘int gsl_bspline_deriv_eval_

--------------------------------------------------------------------------------

## hartreeFock (main program)

 * Solves relativistic Hartree-Fock potential for core + valence states
 * Input taken from a plain text file.
 * See "hartreeFock.in" for minimal input example.
  May re-name this file (e.g., to "filename.txt"), then run as:
    * _$ ./hartreeFock filename.txt_
    * (Otherwise, program will assume file name is 'hartreeFock.in')
 * see doc/ folder for a full list of input options + descriptions


## Documentation

 * Documentation is in doc/ directory (best viewed with a markdown reader or on GitHub). Contains three documents:
 * 01-hartreeFock_input -- How to use the code (input options + descriptions)
 * 02-diracSCAS_method  -- What the code does (description of physics)
 * 03-diracSCAS_code    -- Documentation for code objects/functions etc. [coming "soon" (not soon)]

--------------------------------------------------------------------------------

## Other programs:

### periodicTable

Command-line periodic table, with electron configurations and nuclear data

 * Compiled using the Makefile (run _$make_, must habe 'make' installed)
 * Alternatively, compile with command:
_$g++ src/Physics/AtomData.cpp src/Physics/NuclearData.cpp src/periodicTable.cpp -o periodicTable -I./src/_
 * No other dependencies

Gives info regarding particular element, including Z, default A, and electron configuration.
Takes input in one line from command line.

Usage: (examples)
 * _$./periodicTable_           Prints periodic table
 * _$./periodicTable Cs_        Info for Cs with default A
 * _$./periodicTable Cs 137_    Info for Cs-137
 * _$./periodicTable Cs all_    Info for all available Cs isotopes
 * Note: numbers come from online database, and have some errors,
so should be checked if needed.

 Or, enter 'c' to print list of physics constants
  * _$./periodicTable c_        Prints values for some handy physical constants

Note: ground-state electron configurations are "guessed", and can sometimes be incorrect.

Nuclear radius data mostly comes from:
 * I. Angeli and K. P. Marinova, At. Data Nucl. Data Tables 99, 69 (2013).
https://doi.org/10.1016/j.adt.2011.12.006

Units:
 * r_rms: root-mean-square radius, in fm.
 * c: half-density radius (assuming Fermi nuclear distro)
 * mu: magnetic moment (in nuclear magnetons)


### dmeXSection

 * Calculates the cross-section and event rates for ionisation of atoms
 by scattering of DM particle.
 * Takes in the output of "atomicKernal" (see doc/ for details)
 * Also calculates "observable" event rates, accounting for detector thresholds
 and resolutions (for DAMA/LIBRA and XENON100-like detectors).
 * For definitions/details, see:
   * B.M. Roberts, V.V. Flambaum
[Phys.Rev.D 100, 063017 (2019)](https://link.aps.org/doi/10.1103/PhysRevD.100.063017 "pay-walled");
[arXiv:1904.07127](https://arxiv.org/abs/1904.07127 "free download").
   * B.M.Roberts, V.A.Dzuba, V.V.Flambaum, M.Pospelov, Y.V.Stadnik,
[Phys.Rev.D 93, 115037 (2016)](https://link.aps.org/doi/10.1103/PhysRevD.93.115037 "pay-walled");
[arXiv:1604.04559](https://arxiv.org/abs/1604.04559 "free download").


### wigner

 * Small routine to calculate 3,6,9-j symbols, and Clebsch Gordon coefficients
 * Either give input via command line directly (quote marks required)
   * e.g., _./wigner '<0.5 -0.5, 0.5 0.5| 1 0>'_
   * or e.g., _./wigner '(0.5 1 0.5, -0.5 0 0.5)'_ etc.
 * Or, give an input file, that contains any number of symbols, all on new line
   * e.g., _./wigner -f myInputFile.in_
   * nb: the '-f' flag can be dropped in the '.in' file extension is used
   * Do not use quote marks in input file. Lines marked '!' or '#' are comments
 * 3j symbols must start with '('; 6,9j with '{', and CG with '<' (this is how code knows which symbol to calculate).
 * but, each number can be separated by any symbol (space, comma etc.)
