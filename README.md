# ampsci
_*Atomic Many-body Perturbation theory in the Screened Coulomb Interaction*_

Solves the Dirac equation for single-valence atomic systems using the Hartree-Fock + correlation potential method (based on Dzuba-Flambaum-Sushkov method).
Fully relativistic, includes electron correlations (all-orders screening and hole-particle interaction, finite-nuclear size, Breit interaction, radiative QED effects, and RPA for matrix elements (with structure radiation/renormalisation).
QED is included via the Ginges-Flambaum radiative potential method.
Can solve for continuum states with high energy.
Calculates ionisation cross sections with high values for energy/momentum transfer. Parallelised using openMP.

Designed to be fast, accurate, and easy to use.
The "modules" system (see documentation) makes it simple to add your own routines to use the atomic wavefunctions to calculate whatever properties you may be interested in.

 * The code is on GitHub: [github.com/benroberts999/ampsci](https://github.com/benroberts999/ampsci)
 * This file is best viewed with a markdown reader (or on GitHub)

  [![][doxygen-badge]][docs-url]
  [![][manual-badge]][man-url]
  [![][tests-badge]][actions-url]
  [![][build-badge]][actions-url]

## Documentation

There are four documentation types provided:
 * All documentation available online: [benroberts999.github.io/ampsci/](https://benroberts999.github.io/ampsci/)
 * Also can be found in doc/ directory


 1. README (this document)
    * Brief overview, including compilation instructions (given below)


 2. Input options -- how to run the code
    * _doc/ampsci_input.md_ -- detailed info on all input options
      * Best viewed with a markdown reader or on GitHub
    * See also: _doc/examples/ampsci.in_ -- an example/template input file
    * In _doc/examples/_ there are several example input files, with the expected output; use these to test if everything is working


 3. Physics documentation: _ampsci.pdf_ -- Description of physics/methods used in the code
    * Includes many references to the works where the methods implemented here were developed.
    * Available online: [benroberts999.github.io/ampsci/ampsci.pdf](https://benroberts999.github.io/ampsci/ampsci.pdf)
    * Latex file provided in doc/tex/ampsci.tex
    * If you have latex installed, you can use Makefile to generate the pdf
      * Run '_$make docs_' -- this will create new pdf file: 'doc/ampsci.pdf'


 4. Code documentation -- details on classes/functions in the code
    * Available online: [benroberts999.github.io/ampsci/](https://benroberts999.github.io/ampsci/)
    * Auto-generated (Doxygen), so be wary of mistakes/typos etc.
    * You can generate this by running '_$make doxy_', which produces html documentation (see doc/html/index.html), and a pdf version (doc/documentation.pdf) -- but it's much easier to view this online


--------------------------------------------------------------------------------

## Compilation:

 * Copy "doc/Makefile.example" from doc/ directory to the working directory, and rename to -> "Makefile"
    * _$cp ./doc/Makefile.example ./Makefile_
 * All programs compiled using the Makefile (run _$make_)
 * The file _Makefile_ has some basic compilation options. It's currently set up to work on most linux systems; you may need to change a few options for others (see below)
 * Tested with g++ and clang++ on linux and mac (requires c++17)
    * Works+tested with g++7 and newer
    * Works+tested with clang++6 and newer

Requires GSL (GNU scientific libraries) https://www.gnu.org/software/gsl/, and LAPACK. These must be installed for the code to run (see below).
 * Requires GSL ver 2.0+ (tested with 2.1, 2.6)

 * Quick start (ubuntu). Full dependencies list, run:
   * _$sudo apt install g++ make liblapack-dev libblas-dev libgsl-dev libomp-dev_


### Compilation: Linux:

  * Instructions for ubuntu; similar commands for other flavours
  * Install make: _$sudo apt-get install make_
  * Install GSL libraries: _$sudo apt-get install libgsl-dev_
  * May also need LAPACK/BLAS libraries: _$sudo apt-get install libatlas-base-dev liblapack-dev libblas-dev_
  * Install the compiler: _$sudo apt-get install g++_ and/or _$sudo apt-get install clang++_
  * Then compile by running _$make_ from the ampsci directory
  * To use with openMP (for parallelisation) with clang++, you might have to also install clangs openmp libraries: _$sudo apt install libomp5_ (and perhaps _$sudo apt install libomp-dev_)


### Compilation: MacOS:

  * On mac: use homebrew to install gsl: _$brew install gsl_
  * (homebrew is a package manager; install from _https://brew.sh/_)
  * Seems to work best with the homebrew version of gcc. Install as: _$brew install gcc_
  * Note: you may need to change the compiler from `g++` to `g++-9` (or similar), or update your environment variables, since calling g++ on mac actually calls clang++ by default
  * You might have to tell the compiler how to link to the GSL library; in Makefile:
    - PathForGSL=/usr/local/opt/gnu-scientific-library
  * Then compile by running _$make_ from the ampsci directory
  * Use openMP for parellelisation when using clang++ on mac:
    * If using g++, should work as per normal
    * To use openMP with clang, seem to require the llvm version
    * _$brew install llvm_
    * Then, in the Makefile, set (exact paths may be different for you):
      - CXX=/usr/local/opt/llvm/bin/clang++
      - ExtraInclude=/usr/local/opt/llvm/include/
      - ExtraLink=/usr/local/opt/llvm/lib/
    * Note: For me, I have to set OPT=-O1 for this to work, otherwise takes infinite time to compile.


### Compilation: Windows:

For windows, the easiest way (for me, anyway) is to use the 'windows subsystem for linux' (requires Windows10). Instructions on installation/use here: https://www.roberts999.com/posts/2018/11/wsl-coding-windows-ubuntu.
Then, the compilation + use can proceed as per Linux above.
In theory, should work with MSVC - but no CMAKE file provided (yet)

### Common Compilation errors:

 * **error: unsupported option -fopenmp**
  - openmp (used for parallelisation) is not working. See above for some possible solutions.
  - Quick fix: change '_UseOpenMP=yes_' to '_UseOpenMP=no_' in Makefile

 * **fatal error: gsl/<...>.h: No such file or directory** (or similar)
 * **gsl** related linking/compilation error:
  - Could not find required GSL libraries. Either they are not installed, or you need to link to them
  - 1) Ensure GSL is installed (see above for instructions)
  - 2) If GSL library is not installed in _/usr/local/_, you have to tell the compiler where to find the GSL files. Do this by setting the _PathForGSL_ option in Makefile. Common examples:
    * _PathForGSL=/opt/gsl/2.1/gnu_ # For UQ's getafix cluster
    * _PathForGSL=/usr/local/opt/gnu-scientific-library_ # For my macbook
    * Note: the exact path may differ for you, depending on where GSL was installed

* **error: too few arguments to function ‘int gsl_bspline_deriv_eval**
 - This is because the code is linking to a very old version of GSL. You might need to update GSL. If you have updated GSL (to at least version 2.0) and still get the message, the code is probably linking against the wrong version of GSL; see above to point the compiler to the correct version

* **other:**
 * Sometimes, the compiler will not be able to find the correct libraries (particular, e.g., on clusters). In this case, there are two options in the Makfefile: **ExtraInclude** and **ExtraLink**
   * These add paths to include the correct directories for both -I "includes" (for compilation), and -L link flags (for linking libraries) in Makefile. These can be a little tricky to get right (don't include the -I or -L)

--------------------------------------------------------------------------------

## ampsci (main program)

Run the program with input options from the command line, e.g.:

 * $ ./ampsci filename
   - Runs ampsci with input option specified in file "filename"
   - This is the main way to run program

 * $ ./ampsci <At> <Core> <Valence>
   - For quick use: simple HF calculation. e.g.,
 * $ ./ampsci Cs
    - Runs ampsci for Cs using Hartree Fock (V^N) approximation
 * $ ./ampsci Cs [Xe] 6sd5d
    - Runs ampsci for Cs using Hartree Fock with Xe-like core and valence
      states up to n=6 for s,p-states and n=5 for d-states
 * $ ./ampsci Cs
    - Runs ampsci for Cs using Hartree Fock (V^N) approximation

 * $ ./ampsci -v
    - Prints version info (same as --version)

 * $ ./ampsci -h
    - Print help info, including input options (same as --help)

Output is printed to screen. It's recommended to forward this to a text file.
The input options and the ampsci version details are also printed, so that the
program output contains all required info to exactly reproduce it.
e.g.,
 * ./ampsci input |tee -a outout
   - Runs ampsci using input options in file "input".
   - Output will both be written to screen, and appended to file "output".

 * Input taken from a plain text file (or via command line).
 * Several example input files are given in: _doc/examples/_, along with their expected output; use these to test if everything is working.
 * See _doc/ampsci_input.md_ for a full list of input options + descriptions
 * See _ampsci.pdf_ for a description of the physics, and for references to the works where the methods implemented here were developed.
   * Available on GitHub: [benroberts999.github.io/ampsci/ampsci.pdf](https://benroberts999.github.io/ampsci/ampsci.pdf)


--------------------------------------------------------------------------------

## Other programs:

### periodicTable

Command-line periodic table, with electron configurations, nuclear data, and some physical constants

 * Compiled using the Makefile (run _$make_, must have 'make' installed)
 * Alternatively, compile with command:
   * _$g++ -std=c++17 src/Physics/AtomData.cpp src/Physics/NuclearData.cpp src/periodicTable.cpp -o periodicTable -I./src/_
 * No other dependencies

Gives info regarding particular element, including Z, default A, and electron configuration.
Takes input in one line from command line.

Usage: (examples)
 * _$./periodicTable_         -- Prints periodic table
 * _$./periodicTable Cs_      -- Info for Cs with default A
 * _$./periodicTable Cs 137_  -- Info for Cs-137
 * _$./periodicTable Cs all_  -- Info for all available Cs isotopes
 * Note: numbers come from online database, and have some errors, so should be checked if needed.


 Or, enter 'c' to print list of physics constants
  * _$./periodicTable c_      -- Prints values for some handy physical constants

Note: ground-state electron configurations are "guessed", and can sometimes be incorrect.

Nuclear data mostly comes from:
 * Radius data: I. Angeli and K. P. Marinova, At. Data Nucl. Data Tables 99, 69 (2013).
[doi:10.1016/j.adt.2011.12.006](https://doi.org/10.1016/j.adt.2011.12.006)
 * Magnetic moments: N. Stone, At. Data Nucl. Data Tables 90, 75 (2005).
 [doi:10.1016/j.adt.2005.04.001](https://doi.org/10.1016/j.adt.2005.04.001)
 * Note: data was scraped from the tables, and contains many transcription errors. Provided for convenience, but it is up to you to double check that the values are correct/up-to-date

Units:
 * r_rms: root-mean-square radius, in fm.
 * c: half-density radius (assuming Fermi nuclear distro, with t=2.3)
 * mu: magnetic moment (in nuclear magnetons)


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


### dmeXSection

  * Calculates the cross-section and event rates for ionisation of atoms
  by scattering of DM particle.
    * B.M. Roberts, V.V. Flambaum
 [Phys.Rev.D 100, 063017 (2019)](https://link.aps.org/doi/10.1103/PhysRevD.100.063017 "pay-walled");
 [arXiv:1904.07127](https://arxiv.org/abs/1904.07127 "free download").
    * B.M.Roberts, V.A.Dzuba, V.V.Flambaum, M.Pospelov, Y.V.Stadnik,
 [Phys.Rev.D 93, 115037 (2016)](https://link.aps.org/doi/10.1103/PhysRevD.93.115037 "pay-walled");
 [arXiv:1604.04559](https://arxiv.org/abs/1604.04559 "free download").
  * Old, and has not been tested for a long time; use with caution
  * see _'doc/dmeXSection_input'_ for details


 [docs-url]: https://benroberts999.github.io/ampsci/
 [man-url]: https://benroberts999.github.io/ampsci/ampsci.pdf
 [actions-url]: https://github.com/benroberts999/ampsci/actions
 [build-badge]: https://github.com/benroberts999/ampsci/workflows/Build/badge.svg
 [tests-badge]: https://github.com/benroberts999/ampsci/workflows/Tests/badge.svg
 [doxygen-badge]: https://img.shields.io/badge/documentation-code%20(html)-blue
 [manual-badge]: https://img.shields.io/badge/documentation-physics%20(pdf)-blue
