# Compilation

\brief Compilation instructions for Linux, Mac, and Windows

[[Home](/README.md)]

* Easiest method is to compile using provided Makefile:
* Copy "doc/examples/Makefile" from doc/ directory to the working directory
  * `cp ./doc/examples/Makefile ./`
* All programs compiled using the Makefile (run `make`)
* The file _Makefile_ has some basic compilation options. It's currently set up to work on most linux systems; you may need to change a few options for others (see below)
* Tested with g++ and clang++ on linux and mac

## Dependencies / Requirements

* c++ compiler that supports c++17
  * Tested with clang version 6 and newer; gcc version 7 and newer
  * Also tested with intel [icc 2021.2.0], though this is tested infrequently
* LAPACK and BLAS libraries [netlib.org/lapack/](http://www.netlib.org/lapack/)
* GSL (GNU scientific libraries) [gnu.org/software/gsl/](https://www.gnu.org/software/gsl/) [version 2.0 or newer*]
  * (it _should_ also work with older versions of GSL, but this is not regularly tested and therefore not guarenteed)
* [optional] GNU Make ([gnu.org/software/make/](https://www.gnu.org/software/make/)) - used to compile code
* [optional] OpenMP ([openmp.org/](https://www.openmp.org/)) - used for parallisation
* [optional] git ([git-scm.com/](https://git-scm.com/)) for version tracking and to keep up-to-date with latest version

## Quick-start

This is for ubuntu/linux - for other systems, see below

* Get the code from [GitHub](https://github.com/benroberts999/ampsci), using git:
  * `sudo apt install git`
  * `git clone git@github.com:benroberts999/ampsci.git`
  * or `git clone https://github.com/benroberts999/ampsci.git`
* Or, direct download (without using git, not redcommended):
  * <https://github.com/benroberts999/ampsci/archive/refs/heads/main.zip>
* Install dependencies
  * `sudo apt install g++ liblapack-dev libblas-dev libgsl-dev make libomp-dev`
* Prepare the Makefile (already setup for ubuntu, you may need minor adjustments, see below)
  * `cp ./doc/examples/Makefile ./`
* Compile ampsci using all default options:  
  * `make`
* Run the first example program
  * `cp ./doc/examples/ampsci.in ./`
  * `./ampsci ampsci.in`

--------------------------------------------------------------------------------

## Compilation: Linux

* Instructions for ubuntu; similar commands for other flavours
* Install make: `sudo apt-get install make`
* Install GSL libraries: `sudo apt-get install libgsl-dev`
* May also need LAPACK/BLAS libraries: `sudo apt-get install libatlas-base-dev liblapack-dev libblas-dev`
* Install the compiler: `sudo apt-get install g++` and/or `sudo apt-get install clang++`
* Then compile by running `make` from the ampsci directory
* To use with openMP (for parallelisation) with clang++, you might have to also install clangs openmp libraries: `sudo apt install libomp5` (and/or perhaps `sudo apt install libomp-dev`)

## Compilation: MacOS

* On mac: use _homebrew_ to install gsl: _brew install gsl_
* _homebrew_ is a package manager; install from [https://brew.sh/](https://brew.sh/)
* Seems to work best with the homebrew version of gcc. Install as: `brew install gcc`
* Note: you may need to change the compiler from `g++` to `g++-9` (or similar), or update your environment variables, since calling g++ on mac actually calls clang++ by default
* You might have to tell the compiler how to link to the GSL library; in Makefile:
  * PathForGSL=/usr/local/opt/gnu-scientific-library
* Then compile by running _make_ from the ampsci directory
* Use openMP for parellelisation when using clang++ on mac:
  * If using g++, should work as per normal
  * To use openMP with clang, seem to require the llvm version
  * _brew install llvm_
  * Then, in the Makefile, set (exact paths may be different for you):
    * CXX=/usr/local/opt/llvm/bin/clang++
    * ExtraInclude=/usr/local/opt/llvm/include/
    * ExtraLink=/usr/local/opt/llvm/lib/
  * This seems fragile

## Compilation: Windows

For windows, the easiest way (for me, anyway) is to use the 'windows subsystem for linux' (requires Windows 10+). Instructions on installation/use here: [docs.microsoft.com/en-us/windows/wsl/install](https://docs.microsoft.com/en-us/windows/wsl/install).
Then, the compilation + use can proceed as per Linux above.

--------------------------------------------------------------------------------

## Common Compilation errors

* **error: unsupported option -fopenmp**

* openmp (used for parallelisation) is not working. See above for some possible solutions.
* Quick fix: change '_UseOpenMP=yes_' to '_UseOpenMP=no_' in Makefile

* **fatal error: gsl/<...>.h: No such file or directory** (or similar)
* **gsl** related linking/compilation error:

* Could not find required GSL libraries. Either they are not installed, or you need to link to them
* 1) Ensure GSL is installed (see above for instructions)
* 2) If GSL library is not installed in _/usr/local/_, you have to tell the compiler where to find the GSL files. Do this by setting the _PathForGSL_ option in Makefile. Common examples:
  * _PathForGSL=/opt/gsl/2.1/gnu_ # For UQ's getafix cluster
  * _PathForGSL=/usr/local/opt/gnu-scientific-library_ # For my macbook
  * Note: the exact path may differ for you, depending on where GSL was installed

* **error: too few arguments to function ‘int gsl_bspline_deriv_eval**

* This is because the code is linking to a very old version of GSL. You might need to update GSL. If you have updated GSL (to at least version 2.0) and still get the message, the code is probably linking against the wrong version of GSL; see above to point the compiler to the correct version

* **other:**
* Sometimes, the compiler will not be able to find the correct libraries (particular, e.g., on clusters). In this case, there are two options in the Makfefile: **ExtraInclude** and **ExtraLink**
  * These add paths to include the correct directories for both -I "includes" (for compilation), and -L link flags (for linking libraries) in Makefile. These can be a little tricky to get right (don't include the -I or -L)
