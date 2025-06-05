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
# Set GSL and LAPACK/BLAS library flags. These are typical, but sometimes differ
LIBS=-lgsl -lgslcblas -llapack -lblas

# Friday:
# LIBS=-lgsl -lgslcblas -lopenblas
# Bunya:
# LIBS=-lgsl -lgslcblas -llapack -lopenblas -lgfortran 

################################################################################
# OpenMP library to use. -fopenmp default for GCC, -fopenmp=libomp for clang
# libomp for clang++X requires package libomp-X-dev (e.g., X=15)
# Ignored if UseOpenMP is not 'yes' or 'true'
OMPLIB=-fopenmp
# OMPLIB=-fopenmp=libomp

################################################################################
# Use these to pass in any other compiler/linker arguments
# (Rarely needed)
CARGS=
LARGS=

################################################################################
## None of the below options should need changing
################################################################################
## Set directories for source files (SD), and output built object files (BD)
SD=./src
BD=./build

## c++ standard. must be at least c++17
CXXSTD=-std=c++17

## Build config + options:
include $(BD)/buildOptions.mk

## Build targets (must update if new programs/files are added):
include $(BD)/buildTargets.mk
