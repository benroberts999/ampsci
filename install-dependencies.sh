#!/bin/bash
set -euo pipefail

echo ''
echo 'This will attempt to install all ampsci dependencies.'
echo 'After this, run configure.sh and make to compile ampsci with default options.'
echo 'You may need to update some options in the Makefile, depending on your system.'
echo ''

# Skip confirmation prompt if --yes flag is passed (useful for CI/scripted installs)
if [[ "${1:-}" != "--yes" && "${1:-}" != "-y" ]]; then
  read -r -p "This script will install programs. Do you wish to continue? [y/N] " response
  if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
    echo ""
  else
    exit
  fi
fi

# Detect OS
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    *)          machine="UNKNOWN:${unameOut}"
esac
echo 'Running install dependencies for: '${machine}

if [ ${machine} == 'Mac' ] && [ "$EUID" -eq 0 ]; then
  echo "Do not run this script with sudo on Mac."
  exit 1
fi

################################################################################

# Check if an apt package is installed; if not, add it to missing_apt[]
check_apt() {
  if dpkg -s "$1" &>/dev/null; then
    version=$(dpkg -s "$1" | grep '^Version' | awk '{print $2}')
    echo "  Found $1: $version"
  else
    echo "  Did not find $1: will install"
    missing_apt+=("$1")
  fi
}

install_linux() {
  if ! command -v apt-get &>/dev/null; then
    echo "Detected non-Ubuntu Linux. Please install the following manually:"
    echo "  make, g++, LAPACK, BLAS, GSL"
    echo "e.g.:"
    echo "  dnf:    sudo dnf install make gcc-c++ lapack-devel blas-devel gsl-devel"
    echo "  pacman: sudo pacman -S make gcc lapack blas gsl"
    exit 0
  fi
  missing_apt=()
  check_apt make
  check_apt g++
  check_apt liblapack-dev
  check_apt libblas-dev
  check_apt libgsl-dev
  if [ ${#missing_apt[@]} -eq 0 ]; then
    echo 'All dependencies already installed.'
  else
    set -x
    sudo apt-get update &&
    sudo apt-get install --yes "${missing_apt[@]}"
  fi
}

################################################################################

install_homebrew() {
  if ! command -v brew &>/dev/null; then
    read -r -p "Homebrew not installed. Would you like to install it? [y/N] " response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
      set -x
      /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    else
      exit
    fi
  fi
}

################################################################################

# Check Mac dependencies directly (without needing brew);
# plain g++ on Mac is Apple Clang, so look for a versioned g++ (same as configure.sh)
check_mac_deps() {
  missing_brew=()
  found_gcc=0
  for ver in $(seq 30 -1 7); do
    if [[ -x "$(command -v "g++-${ver}")" ]]; then
      echo "  Found gcc: g++-${ver}"
      found_gcc=1
      break
    fi
  done
  if (( found_gcc == 0 )); then
    echo "  Did not find gcc: will install"
    missing_brew+=(gcc)
  fi
  # gsl-config is installed alongside GSL and is the standard way to detect it
  if [[ -x "$(command -v gsl-config)" ]]; then
    echo "  Found gsl: $(gsl-config --version)"
  else
    echo "  Did not find gsl: will install"
    missing_brew+=(gsl)
  fi
}

install_mac() {
  check_mac_deps
  if [ ${#missing_brew[@]} -eq 0 ]; then
    echo 'All dependencies already installed.'
  else
    install_homebrew
    set -x
    brew install "${missing_brew[@]}"
  fi
}

################################################################################

if [ ${machine} == 'Linux' ]; then
  install_linux
elif [ ${machine} == 'Mac' ]; then
  install_mac
else
  echo 'Unsupported system for setup - requires manual setup'
  exit
fi
