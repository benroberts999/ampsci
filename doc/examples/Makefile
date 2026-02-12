################################################################################
# Current options are set up for standard linux
# You may need to tweak the options
################################################################################

## Which compiler: (g++, clang++)
CXX ?= g++

## Build mode (changes warnings + optimisation level): release/dev/debug
MODE ?= release

## Speed up compilation by doing in parallel
DEFAULT_MFLAGS = -j6

## Optimisation flags
OPT ?= -O3

################################################################################
## Path to the GSL library. Usually, this can be left blank. 
## Exact path will depend on where GSL has been installed. 
GSL_PATH ?=

## Other common place:
# GSL_PATH ?= /opt/gsl/2.1/gnu
## Intel Mac: 
# GSL_PATH ?= /usr/local/opt/gnu-scientific-library
## M1/M2 Mac: 
# GSL_PATH ?= /opt/homebrew/Cellar/gsl/2.7.1/

################################################################################
## Set GSL and LAPACK/BLAS library flags. These are typical, but sometimes differ
LIBS ?= -lgsl -lgslcblas -llapack -lblas

## Friday:
# LIBS ?= -lgsl -lgslcblas -lopenblas

## Bunya:
# LIBS ?= -lgsl -lgslcblas -llapack -lopenblas -lgfortran 

################################################################################
## OpenMP library to use. Comment out or leave blank for no OpenMP
OMPLIB ?= -fopenmp

## Usually '-fopenmp' default for GCC, -fopenmp=libomp for clang
# OMPLIB ?= -fopenmp=libomp

################################################################################
# Use these to pass in any other compiler/linker arguments (Rarely needed)
CARGS ?=
LARGS ?=

################################################################################
## The following must be inlcuded
## Build config + options:
include ./src/buildOptions.mk
include ./src/buildRules.mk
include ./src/buildTargets.mk
