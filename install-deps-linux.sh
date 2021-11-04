#!/bin/bash
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
apt install --yes libomp-dev
echo ''
echo 'To compile AMPSCI, copy the example makefile to the working directory:'
echo '--  $ cp ./doc/Makefile.example ./Makefile'
echo 'Then, run make:'
echo '--  $ make'
