\page compilation Compilation Details
\ingroup getting_started
\brief Detailed compilation instructions for Linux, Mac, and Windows

## Dependencies / Requirements

* C++ compiler that supports C++17
  * Tested with clang version 6 and newer; gcc version 7 and newer
  * Tested with g++ and clang++ on linux and mac
  * Also tested with Intel [icc 2021.9.0], though this is tested infrequently
* LAPACK and BLAS libraries [netlib.org/lapack/](http://www.netlib.org/lapack/)
* GSL (GNU scientific libraries) [gnu.org/software/gsl/](https://www.gnu.org/software/gsl/) [version 2.0 or newer*]
  * (it _should_ also work with older versions of GSL, but this is not regularly tested and therefore not guaranteed)
  * Updated to work with newer GSL version 2.8; still works with older versions too
* [optional] GNU Make ([gnu.org/software/make/](https://www.gnu.org/software/make/)) - used to compile code
* [optional] OpenMP ([openmp.org/](https://www.openmp.org/)) - used for parallelisation
* [optional] git ([git-scm.com/](https://git-scm.com/)) for version tracking and to keep up-to-date with latest version
* The shell script `install-dependencies.sh` will attempt to automatically install the required dependencies

--------------------------------------------------------------------------------

## Compilation: Linux

* Instructions for ubuntu; similar commands for other flavours
* Install make: `sudo apt-get install make`
* Install GSL libraries: `sudo apt-get install libgsl-dev`
* May also need LAPACK/BLAS libraries: `sudo apt-get install libatlas-base-dev liblapack-dev libblas-dev`
* Install the compiler: `sudo apt-get install g++` and/or `sudo apt-get install clang++`
* Then compile by running `make` from the ampsci directory
* To use with openMP (for parallelisation) with clang++, you might have to also install clang openmp libraries: `sudo apt install libomp-dev`
  * You'll often need the exact versions to match. So, if using `clang++-15`, install `libomp-dev-15`.

## Compilation: MacOS (intel chip)

* On mac: use _homebrew_ to install gsl: _brew install gsl_
* _homebrew_ is a package manager; install from [https://brew.sh/](https://brew.sh/)
* Seems to work best with the homebrew version of gcc. Install as: `brew install gcc`
* Note: you may need to change the compiler from `g++` to `g++-9` (or similar), or update your environment variables, since calling g++ on mac actually calls clang++ by default
* You might have to tell the compiler how to link to the GSL library; in Makefile:
  * `GSL_PATH=/usr/local/opt/gnu-scientific-library`
* Then compile by running _make_ from the ampsci directory
* Use openMP for parellelisation when using clang++ on mac:
  * If using g++, should work as per normal
  * To use openMP with clang, seem to require the llvm version
  * _brew install llvm_
  * Then, in the Makefile, set (exact paths may be different for you):
    * `CXX = /usr/local/opt/llvm/bin/clang++`
    * `CARGS += -I/usr/local/opt/llvm/include/`
    * `LARGS += -l/usr/local/opt/llvm/lib/`
  * This seems fragile

## Compilation: MacOS (M1/M2/apple silicon chip)

* Mostly the same as above, but some of the libraries are installed to different directories by default. In particular:
* `GSL_PATH=/opt/homebrew/Cellar/gsl/2.7`

## Compilation: Windows

* For windows, the easiest way (for me, anyway) is to use the 'windows subsystem for linux' (requires Windows 10+).
* Instructions on installation/use here: [docs.microsoft.com/en-us/windows/wsl/install](https://docs.microsoft.com/en-us/windows/wsl/install).
* Then, the compilation + use can proceed as per Linux above.
* You may have some luck getting the MSVC compiler working, but linking to libraries will be difficult

--------------------------------------------------------------------------------

## Compilation Options

Various basic options can be specified in the `Makefile`:

* Set the compiler you want to use:

```make
## Which compiler: (g++, clang++)
CXX = g++
```

* Build modes:
  * Choose one of `release/dev/debug`
  * `release` (default) -- turns off warnings and checks. This is fine for main branch if you haven't modified anything, but should not be used if you are modiying the code
  * `dev` -- please use `dev` mode, which turns on checks and warning, if you are modyfying the code. This will help to catch any errors or issues
  * `debug` -- will be slow, but turns on some debugging options to make debugging easier

```make
## Build mode (changes warnings + optimisation level): release/dev/debug
MODE ?= release
```

* Path to the GSL libraries. This can be left blank the majority of the time. On some systems if GSL is not found natively, or if the incorrect labrary is being found, you might need to specify it
  * e.g., on newer M1/M2 macs, you may need `/opt/homebrew/Cellar/gsl/2.7.1/`
  * intel mac: `/usr/local/opt/gnu-scientific-library`
  * the `configure.sh` script should have done this for you

```make
## Path to the GSL library. Usually, this can be left blank. 
## Exact path will depend on where GSL has been installed. 
GSL_PATH ?=
```

* Libraries. Again, usually the default is fine.
  * On some systems, different libraries are used (e.g., openblas vs. blas)
  * In particular, on Firday cluster: `-lgsl -lgslcblas -lopenblas`
  * On Bunya: `-lgsl -lgslcblas -llapack -lopenblas -lgfortran`

```make
## Set GSL and LAPACK/BLAS library flags. These are typical, but sometimes differ
LDLIBS ?= -lgsl -lgslcblas -llapack -lblas
```

* OpenMP library
* Comment out or leave option blank to turn off OpenMP multithreading

```make
## OpenMP library to use. Comment out or leave blank for no OpenMP
OMPLIB ?= -fopenmp
```

* Extra compiler/linker arguments: usually blank, but you may add to these, which are passed as compiler and linker arguments

```make
# Use these to pass in any other compiler/linker arguments (Rarely needed)
CARGS ?=
LARGS ?=
```
