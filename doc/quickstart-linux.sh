#!/bin/bash
# Check if run as root:
if [ "$EUID" -ne 0 ]
  then echo "Please run using sudo"
  exit
fi
# Install all the dependencies for AMPSCI
echo ''
echo 'Install all basic dependencies for AMPSCI'
echo ''
echo 'Install make:'
apt install --yes make &&
echo ''
echo 'Install c++ compiler (g++) - you may use another c++ compiler:'
apt install --yes g++ &&
echo ''
echo 'Install LAPACK (Linear Algebra package):'
apt install --yes liblapack-dev &&
echo ''
echo 'Install BLAS (Basic Linear Algebra Subprograms):'
apt install --yes libblas-dev &&
echo ''
echo 'Install GSL (GNU Scientific Libraries):'
apt install --yes libgsl-dev &&
echo ''
echo 'Install OpenMP libraries (optional):'
apt install --yes libomp-dev  &&
echo ''
echo 'Compile AMPSCI: copy example makefile to the working directory:'
echo 'Then, run make:'
cp ./doc/examples/Makefile ./  &&
cp ./doc/examples/ampsci.in ./ampsci.in &&
make