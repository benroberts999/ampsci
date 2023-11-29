#!/bin/bash

echo ''
echo 'This attempt to compile ampsci using'
echo 'the default options. You may need to update some options in the Makefile,'
echo 'depending on your system.'
echo 'You may need to run install-dependencies.sh first, to install dependencies'
echo ''

if [[ "$1" != "--yes" ]]; then
  read -r -p "This script will compile programs. Do you wish to continue? [y/N] " response
  if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
  then
    echo ""
  else
    exit
  fi
else
  echo ""
fi

# Figure out which operating system is used:
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    *)          machine="UNKNOWN:${unameOut}"
esac
echo 'Running setup for: '${machine}
if [ ${machine} == "UNKNOWN:${unameOut}" ]; then
  echo 'Unsupported system for setup - requires manual setup'
  exit
fi

# Copy across example Makefile and input file (if not already there)
cp -n ./doc/examples/Makefile ./
cp -n ./doc/examples/ampsci.in ./ampsci.in

# Locate GCC compiler: find newest installed first:
gccver=14
while [[ $gccver -gt 6 ]]; do
  if [ -x "$(command -v "g++-"${gccver})" ]; then
    echo "Compiler Found:" $(command which "g++-"${gccver})
    cxx=$(command which "g++-"${gccver})
    sed -i '' -e "s@CXX=.*@CXX=${cxx}@g" Makefile 2> /dev/null
    break
  fi
  ((gccver--))
done
# If didn't find specific version, try g++
if ! [ $gccver -gt 6 ]; then
  if [ -x "$(command -v "g++")" ]; then
    echo "Compiler Found:" "$(command -v "g++")"
    cxx="g++"
    sed -i '' -e "s@CXX=.*@CXX=${cxx}@g" Makefile 2> /dev/null
  else
    echo "Compatible gcc compiler not found."
    echo "Run: install-dependencies.sh, or"
    if [ ${machine} == 'Linux' ]; then
      echo "Run: sudo apt install g++"
    elif [ ${machine} == 'Mac' ]; then
      echo "Run: brew install gcc"
    fi
    echo "Or proceed with manual compilation for other compilers"
    exit
  fi
fi

# Locate GSL libraries:
if [ -x "$(command -v "gsl-config")" ]; then
  gslver=$(gsl-config --version)
  gslpath=$(gsl-config --prefix)
  echo "GSL found version:" $gslver
  echo "GSL Libraries found:" $gslpath
  sed -i '' -e "s@PathForGSL=.*@PathForGSL=${gslpath}@g" Makefile 2> /dev/null
else
  echo "GSL libraries not installed"
  echo "Run: install-dependencies.sh, or"
  if [ ${machine} == 'Linux' ]; then
    echo "Run: sudo apt install libgsl-dev"
  elif [ ${machine} == 'Mac' ]; then
    echo "Run: brew install gsl"
  fi
  exit
fi

# Attempt to compile ampsci
echo ''
echo "Attempt to compile ampsci:"
echo ''

make
if [ $? -eq 0 ]; then
  echo ''
  echo "ampsci succesfully compiled"
else
  echo ''
  echo "ampsci compilation failed"
  echo "Check: "
  echo "  - LAPACK/BLAS installed and linking"
  echo "  - blas libraries (blas vs openblas), openmp (set to false)"
  echo "  - openmp (set UseOpenMP=no if not installed)"
  echo ""
fi
