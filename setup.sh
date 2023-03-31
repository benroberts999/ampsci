#!/bin/bash

echo 'This will install all dependencies, and attempt to compile ampsci using'
echo 'the default options. You may need to update some options in the Makefile,'
echo 'depending on your system'

if [[ "$1" != "--yes" ]]; then
  read -r -p "This script will install programs. Do you wish to continue? [y/N] " response
  if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
  then
    echo ""
  else
    exit
  fi
else
  echo ""
fi

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    *)          machine="UNKNOWN:${unameOut}"
esac
echo 'Running setup for: '${machine}

cp -n ./doc/examples/Makefile ./

cp -n ./doc/examples/ampsci.in ./ampsci.in

if [ ${machine} == 'Linux' ]; then
  # Check if run as root:
  if [ "$EUID" -ne 0 ]; then
    echo "Please run using sudo"
    exit
  fi
  set -x
  apt-get update &&
  apt-get install --yes make g++ liblapack-dev libblas-dev libgsl-dev libomp-dev &&
  make
elif [ ${machine} == 'Mac' ]; then
  # check if homebrew is installed:
  which -s brew
  if [[ $? != 0 ]] ; then
        # If not installed, ask if user would like to install it:
        read -r -p "Homebrew not installed Would you like to install it? [y/N] " response
        if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
          # Install Homebrew
          set -x
          /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        else
          exit
        fi
  else
    brew update
  fi
  set -x
  brew install gcc@12 &&
  brew install gsl &&
  # nb: path to gsl is different for new M1 mac. 
  # I don't have an M1 mac, so this isn't tested
  if [[ $(arch) == 'arm64' ]]; then
    echo "M1 silicon mac"
    make CXX=g++-12 PathForGSL=/opt/homebrew/Cellar/gsl/2.7
  else
    echo "Intel silicon mac"
    make CXX=g++-12 PathForGSL=/usr/local/opt/gsl
    # nb: simply using g++ here will link to clang instead
  fi
else
  echo 'Unsupported system for setup - requires manual setup'
  exit
fi


