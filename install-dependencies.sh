#!/bin/bash

echo ''
echo 'This will attempt to install all ampsci dependencies.'
echo 'After this, run setup.sh to compile ampsci with default options.'
echo 'You may need to update some options in the Makefile, depending on your system.'
echo ''

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
echo 'Running install dependencies for: '${machine}


if [ ${machine} == 'Linux' ]; then
  # Check if run as sudo:
  if [ "$EUID" -ne 0 ]; then
    echo "Please run using sudo"
    exit
  fi
  set -x
  apt-get update &&
  apt-get install --yes make g++ liblapack-dev libblas-dev libgsl-dev libomp-dev
elif [ ${machine} == 'Mac' ]; then
  # check if homebrew is installed:
  which -s brew
  if [[ $? != 0 ]] ; then
        # If not installed, ask if user would like to install it:
        read -r -p "Homebrew not installed. Would you like to install it? [y/N] " response
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
  brew install gcc &&
  brew install gsl
else
  echo 'Unsupported system for setup - requires manual setup'
  exit
fi


