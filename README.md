# diracSCAS: Relativistic, self-consistent atomic structure code.

Solves the Dirac equation for single-valence atomic systems using the Hartree-Fock + second-order correlation potential method.
Fully relativistic, includes finite-nuclear size, Breit interaction, QED effects, RPA for matrix elements, and can solve for continuum states (energy normalisation).
QED corrections included via the Ginges-Flambaum radiative potential method.
Calculates ionisation cross sections with high values for energy/momentum transfer.

 * The code is on GitHub: [github.com/benroberts999/diracSCAS](https://github.com/benroberts999/diracSCAS)

 * This file is best viewed with a markdown reader (or on GitHub)


## Documentation

 * Easiest is to view online: [benroberts999.github.io/diracSCAS/](https://benroberts999.github.io/diracSCAS/)
 * Also available in doc/

There are three documentation types provided:

 1. Input options -- how to run the code
    * _doc/diracSCAS_input.md_ -- detailed info on all input options
      * Best viewed with a markdown reader or on GitHub
    * See also: _doc/diracSCAS.in.example_ -- an example input file for Cs
      * copy to main directory + remove the '.example'
      * _$cp ./doc/diracSCAS.in.example ./diracSCAS.in_


 2. Physics documentation: _diracSCAS.pdf_ -- Description of physics/methods used in the code
    * Includes many references to the works where the methods implemented here were developed.
    * Available online: [benroberts999.github.io/diracSCAS/diracSCAS.pdf](https://benroberts999.github.io/diracSCAS/diracSCAS.pdf)
    * Latex file provided in doc/tex/diracSCAS.tex
    * If you have latex installed, you can use Makefile to generate the pdf
      * Run '_$make docs_' -- this will create new pdf file: 'doc/diracSCAS.pdf'


 3. Code documentation -- details on classes/function in the code
    * Available online: [benroberts999.github.io/diracSCAS/](https://benroberts999.github.io/diracSCAS/)
    * Auto-generated (Doxygen), so be wary of mistakes/typos etc.
    * You can generate this by running '_$make doxy_', which produces html documentation (see doc/html/index.html), and a pdf version (doc/documentation.pdf) -- but it's much easier to view this online


--------------------------------------------------------------------------------

## Compilation:

 * Copy "doc/Makefile.example" from doc/ directory to the working directory, and rename to -> "Makefile"
    * _$cp ./doc/Makefile.example ./Makefile_
 * All programs compiled using the Makefile (run _$make_)
 * The file _Makefile_ has some basic compilation options. It's currently set up to work on most linux systems; you may need to change a few options for others (see below)
 * Tested with g++ and clang++ on linux and mac (requires c++17)

Note: makes use of GSL libraries (requires ver 2.0+, tested with 2.1, 2.6) https://www.gnu.org/software/gsl/, and LAPACK. These must be installed for the code to run (see below).


### Compilation: Linux:

  * Instructions for ubuntu; similar commands for other flavours
  * Install make: _$sudo apt-get install make_
  * Install GSL libraries: _$sudo apt-get install libgsl-dev_
  * May also need LAPACK/BLAS libraries: _$sudo apt-get install libatlas-base-dev liblapack-dev libblas-dev_
  * Install the compiler: _$sudo apt-get install g++_ and/or _$sudo apt-get install clang++_
  * Then compile by running _$make_ from the diracSCAS directory
  * To use with openMP (for parallelisation) with clang++, you might have to also install clangs openmp libraries: _$sudo apt install libomp5_ (and perhaps _$sudo apt install libomp-dev_)


### Compilation: MacOS:

  * On mac: use homebrew to install gsl: _$brew install gsl_
  * (homebrew is a package manager; install from _https://brew.sh/_)
  * Seems to work best with the homebrew version of gcc. Install as: _$brew install gcc_
  * Note: you may need to change the compiler from `g++` to `g++-9` (or similar), or update your environment variables, since calling g++ on mac actually calls clang++ by default
  * You might have to tell the compiler how to link to the GSL library; see below
  * Then compile by running _$make_ from the diracSCAS directory


### Compilation: Windows:

For windows, the easiest way (for me, anyway) is to use the 'windows subsystem for linux' (requires Windows10). Instructions on installation/use here: https://www.roberts999.com/posts/2018/11/wsl-coding-windows-ubuntu.
Then, the compilation + use can proceed as per Linux above.

### Compilation: Other:

 * NOTE: If you get the following error message on compile, change '_UseOpenMP=yes_' to '_UseOpenMP=no_' in Makefile:
   * **error: unsupported option -fopenmp**

 * **Linking to GSL**: If GSL library is not installed in _/usr/local/_, you have to tell the compiler where to find the GSL files. Do this by setting the _PathForGSL_ option in Makefile. Common examples:
   * _PathForGSL=/opt/gsl/2.1/gnu_ # For UQ's getafix cluster
   * _PathForGSL=/usr/local/opt/gnu-scientific-library_ # For my macbook
   * Note: the exact path may differ for you, depending on where GSL was installed

 * Sometimes, the compiler will not be able to find the correct libraries (particular, e.g., on clusters). In this case, there are two options in the Makfefile: **ExtraInclude** and **ExtraLink**
   * These add paths to include the correct directories for both -I "includes" (for compilation), and -L link flags (for linking libraries) in Makefile. These can be a little tricky to get right (don't include the -I or -L)

 * NOTE: If you get the following error, it is because the code is linking to a very old version of GSL. You might need to update GSL. If you have updated GSL (to at least version 2.0) and still get the message, the code is probably linking against the wrong version of GSL; see above to point the compiler to the correct version
   * **error: too few arguments to function ‘int gsl_bspline_deriv_eval**

--------------------------------------------------------------------------------

## diracSCAS (main program)

 * Solves the Dirac equation for atomic systems using the Hartree-Fock + second-order correlation potential method.
 * Input taken from a plain text file.
 * An example input file is included: doc/diracSCAS.in.example
    * e.g.: _$ cp ./doc/diracSCAS.in.example ./diracSCAS.in_
 * You may re-name this file (e.g., to "filename.txt"), then run as:
    * _$ ./diracSCAS filename.txt_
    * If no input filename is given, program will assume input filename is 'diracSCAS.in':
 * Note: input file uses c++-like format + line comments; tell your editor that the file is a cpp file to get nice colourisation and auto commenting
 * See _doc/diracSCAS_input.md_ for a full list of input options + descriptions
 * See _diracSCAS.pdf_ for a description of the physics, and for references to the works where the methods implemented here were developed.
   * Available on GitHub: [benroberts999.github.io/diracSCAS/diracSCAS.pdf](https://benroberts999.github.io/diracSCAS/diracSCAS.pdf)

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
   * B.M. Roberts, V.V. Flambaum
[Phys.Rev.D 100, 063017 (2019)](https://link.aps.org/doi/10.1103/PhysRevD.100.063017 "pay-walled");
[arXiv:1904.07127](https://arxiv.org/abs/1904.07127 "free download").
   * B.M.Roberts, V.A.Dzuba, V.V.Flambaum, M.Pospelov, Y.V.Stadnik,
[Phys.Rev.D 93, 115037 (2016)](https://link.aps.org/doi/10.1103/PhysRevD.93.115037 "pay-walled");
[arXiv:1604.04559](https://arxiv.org/abs/1604.04559 "free download").
 * see _'doc/dmeXSection_input'_ for details


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
