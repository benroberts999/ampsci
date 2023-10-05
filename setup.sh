#!/bin/bash

echo 'This attempt to compile ampsci using'
echo 'the default options. You may need to update some options in the Makefile,'
echo 'depending on your system.'
echo 'You may need to run install-dependencies.sh first, to install dependencies'

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
  make
elif [ ${machine} == 'Mac' ]; then
  set -x
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


