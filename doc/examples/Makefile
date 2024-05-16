################################################################################
# Copy this file into the working directory, and rename to "Makefile"
# Then, compile the code simply by typing 'make' from command line
# Current options are set up for linux; you may need to tweak the options
################################################################################

## Which compiler: (g++, clang++) [no spaces]
CXX=g++

## Use OpenMP (parallelisation) yes/no:
UseOpenMP=yes

## Build mode (changes warnings + optimisation level): release/dev/debug
Build=release

## Compile in parallel with n cores (Faster, but makes reading errors difficult)
ParallelBuild=6

## Optional: set directory for executables (by default: current directory)
XD=.

################################################################################
## Path to the GSL library. Try first with this blank. Exact path will depend
## on where GSL has been installed. Usually, this can be left blank.
## nb: do not leave a trailing space after the line
## Some common places include '/opt/gsl/2.1/gnu'
## Intel Mac: '/usr/local/opt/gnu-scientific-library', 
## M1/M2 mac: '/opt/homebrew/Cellar/gsl/2.7.1/'

PathForGSL=

################################################################################
# By default, uses -lblas (set =blas). For openblas, set =openblas.
# On some systems, openblas is available
BLASLIB=blas

################################################################################
# OpenMP library to use. -fopenmp default for GCC, -fopenmp=libomp for clang
# libomp for clang++X requires package libomp-X-dev (e.g., X=15)
OMPLIB=-fopenmp
# OMPLIB=-fopenmp=libomp

################################################################################
## If compiler cannot find correct libraries/headers, add the paths here.
## (Adds to -I and -L on compile and link; don't include the "-L" or "-I" here)
## Usually, these will be blank.
ExtraInclude=
ExtraLink=
ExtraFlags=
## Required on bunya:
#ExtraFlags=-lgfortran 

# To use openmp with the llvm version clang++ on mac:
#ExtraInclude=/usr/local/opt/llvm/include/
#ExtraLink=/usr/local/opt/llvm/lib/

################################################################################
# Use these to pass in any other compiler/linker arguments
# Any other arguments
CARGS=
LARGS=

################################################################################
## None of the below options should need changing
################################################################################
## Set directories for source files (SD), and output object files (OD)
SD=./src
BD=./build

## c++ standard. must be at least c++17
CXXSTD=-std=c++17

## Build config + options:
include $(BD)/buildOptions.mk

## Build targets (must update if new programs/files are added):
include $(BD)/buildTargets.mk
