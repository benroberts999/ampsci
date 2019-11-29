## Which compiler: (g++, clang++) [no spaces]
CXX=g++
#CXX=clang++

## Use OpenMP (parallelisation) yes/no:
UseOpenMP=yes
#UseOpenMP=no

## Build mode (changes warnings + optimisation level)
#Build=release
Build=dev
#Build=debug

## Optional: set directory for executables (by default: current directory)
XD=.

################################################################################
## None of the below options should need changing
################################################################################
## Set directories for source files (SD), and output object files (OD)
SD=./src
BD=./build

## c++ standard. must be at least c++14 ==> 17
CXXSTD=-std=c++17
#CXXSTD=-std=c++14

## Build config + options:
include $(BD)/buildOptions.mk

## Build targets (must update if new programs/files are added):
include $(BD)/buildTargets.mk
