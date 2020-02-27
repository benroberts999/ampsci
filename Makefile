## Which compiler: (g++, clang++) [no spaces]
CXX=g++
#CXX=clang++

## Use OpenMP (parallelisation) yes/no:
UseOpenMP=yes
#UseOpenMP=no

## Build mode (changes warnings + optimisation level)
Build=release
#Build=dev
#Build=debug

## Compile in parallel with n cores (Faster, but makes reading errors difficult)
ParallelBuild=6

## Optional: set directory for executables (by default: current directory)
XD=.

################################################################################
# If compiler cannot find correct libraries/headers, add the paths here.
# (Adds to -I and -L on compile and link; don't include the "-L" or "-I" here)
# Current defaults are for UQ's getafix to find correct GSL headers/libraries

#ExtraInclude=/opt/gsl/2.1/gnu/include
#ExtraLink=/opt/gsl/2.1/gnu/lib/
#ExtraFlags=

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
