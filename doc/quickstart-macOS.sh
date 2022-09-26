#!/bin/bash
# Install all the dependencies for AMPSCI
echo ''
echo 'Install all basic dependencies for AMPSCI'
echo ''
echo 'Install Homebrew:'
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" &&
brew install gcc@11 &&
brew install gsl &&
echo ''
echo 'Compile AMPSCI: copy example makefile to the working directory:'
echo 'Then, run make:'
cp ./doc/examples/Makefile ./  &&
cp ./doc/examples/ampsci.in ./ampsci.in &&
make make CXX=g++-11 UseOpenMP=yes PathForGSL=/usr/local/opt/gnu-scientific-library