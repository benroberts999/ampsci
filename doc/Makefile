################################################################################
# Example makefile for ampsci. Copy to working directory.
# Currently set up for standard linux; you may need to tweak options.
# The configure.sh script will (attempt to) do this for you.
################################################################################

## Which compiler: (g++, clang++)
CXX = g++

## Build mode (changes warnings + optimisation level): release/dev/debug
## Change to dev if you plan to edit the code
MODE ?= release

## Speed up compilation by doing in parallel
DEFAULT_MFLAGS = -j6

## Optimisation flags
OPT ?= -O3

################################################################################
## Path to the GSL library. Usually, this can be left blank. 
## Exact path will depend on where GSL has been installed. 
GSL_PATH ?=

## Ubuntu/Linux: usually fine to leave blank
## Mac: use `brew --prefix gsl` to find path

################################################################################
## Set GSL and LAPACK/BLAS library flags. These are typical, but sometimes differ
LDLIBS ?= -lgsl -lgslcblas -llapack -lblas

## Friday:
# LDLIBS ?= -lgsl -lgslcblas -lopenblas

## Bunya:
# LDLIBS ?= -lgsl -lgslcblas -llapack -lopenblas -lgfortran 

################################################################################
## OpenMP library to use. Comment out or leave blank for no OpenMP
OMPLIB ?= -fopenmp

################################################################################
# Use these to pass in any other compiler/linker arguments (Rarely needed)
CXXFLAGS +=
LDFLAGS +=

################################################################################
## The following must be inlcuded
## Build config + options:
include ./src/buildOptions.mk
include ./src/buildRules.mk
include ./src/buildTargets.mk
