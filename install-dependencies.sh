#!/bin/bash
set -euo pipefail

# Parse flags: --yes/-y to skip confirmation, --test to simulate all deps missing
AUTO_YES=false
TEST_MODE=false
for arg in "$@"; do
  case "$arg" in
    --yes|-y) AUTO_YES=true ;;
    --test)   TEST_MODE=true ;;
  esac
done

echo ''
echo 'Checking ampsci dependencies...'
if $TEST_MODE; then echo '  (test mode: simulating all dependencies missing)'; fi
echo ''

# Detect OS
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    *)          machine="UNKNOWN:${unameOut}"
esac
echo "Platform: ${machine}"

# Running as root is unnecessary and risky -- the script calls sudo itself when needed
if [ "$EUID" -eq 0 ]; then
  echo "You should not run this script with sudo. We will ask nicely when we need it."
  exit 1
fi

# Minimum required compiler major versions
MIN_GXX_VERSION=7
MIN_CLANGXX_VERSION=6

# Returns the major version number of a compiler (e.g. "11" from "g++ 11.4.0")
compiler_major_version() {
  "$1" -dumpversion 2>/dev/null | cut -d. -f1
}

################################################################################
# apt (Ubuntu/Debian)
################################################################################

# Check if an apt package is installed; if not, report and add to missing_apt[]
check_apt() {
  if ! $TEST_MODE && dpkg -s "$1" &>/dev/null; then
    version=$(dpkg -s "$1" | grep '^Version' | awk '{print $2}')
    echo "  Found $1: $version"
  else
    echo "  Missing: $1"
    missing_apt+=("$1")
  fi
}

# Check for a suitable C++ compiler and append to the named array if not found.
# Takes the array name as a string and uses eval to append -- avoids declare -n
# so this works on bash 3 (macOS default) and old Linux systems.
check_compiler() {
  local missing_ref=$1
  if ! $TEST_MODE && command -v g++ &>/dev/null; then
    ver=$(compiler_major_version g++)
    if (( ver >= MIN_GXX_VERSION )); then
      echo "  Found g++: $(g++ --version | head -1)"
      echo "  -- g++ version: ${ver} (OK)"
    else
      echo "  Found g++ v${ver}, but need >= ${MIN_GXX_VERSION}: will install"
      eval "${missing_ref}+=(g++)"
    fi
  elif ! $TEST_MODE && command -v clang++ &>/dev/null; then
    ver=$(compiler_major_version clang++)
    if (( ver >= MIN_CLANGXX_VERSION )); then
      echo "  Found clang++: $(clang++ --version | head -1) -- skipping g++"
      echo "  -- clang++ version: ${ver} (OK)"
    else
      echo "  Found clang++ v${ver}, but need >= ${MIN_CLANGXX_VERSION}: will install g++"
      eval "${missing_ref}+=(g++)"
    fi
  else
    echo "  Missing: g++"
    eval "${missing_ref}+=(g++)"
  fi
}

# Check deps using only portable tools (command -v, gsl-config, ldconfig).
# Used on non-apt Linux where dpkg is not available -- reports but cannot install.
check_linux_deps_generic() {
  if command -v make &>/dev/null; then
    echo "  Found make: $(make --version | head -1)"
  else
    echo "  Missing: make"
  fi
  # __dummy__ is a dummy array name -- we only want the printed output here
  check_compiler __dummy__
  # gsl-config is installed alongside GSL and is the standard way to detect it
  if ! $TEST_MODE && command -v gsl-config &>/dev/null; then
    echo "  Found gsl: $(gsl-config --version)"
  else
    echo "  Missing: gsl"
  fi
  # ldconfig -p lists all cached shared libraries -- best portable check for LAPACK/BLAS
  for lib in liblapack libblas; do
    if ! $TEST_MODE && ldconfig -p 2>/dev/null | grep -q "${lib}"; then
      echo "  Found: ${lib}"
    else
      echo "  Missing (or could not verify): ${lib}"
    fi
  done
}

# Check deps using apt/dpkg (Ubuntu/Debian only)
check_linux_deps() {
  missing_apt=()
  check_apt make
  check_compiler missing_apt
  check_apt liblapack-dev
  check_apt libblas-dev
  check_apt libgsl-dev
}

install_linux() {
  if ! command -v apt-get &>/dev/null; then
    # Not an apt-based system -- check what we can, then give manual instructions
    echo "Non-apt Linux detected -- checking dependencies anyway:"
    check_linux_deps_generic
    echo ''
    echo "Cannot install automatically. Install missing dependencies manually, e.g.:"
    echo "  dnf:    sudo dnf install make gcc-c++ lapack-devel blas-devel gsl-devel"
    echo "  pacman: sudo pacman -S make gcc lapack blas gsl"
    exit 0
  fi
  check_linux_deps
  if [ ${#missing_apt[@]} -eq 0 ]; then
    echo ''
    echo 'All dependencies already installed.'
    exit 0
  fi
  echo ''
  echo "Need to install: ${missing_apt[*]}"
  echo "  Command: sudo apt-get install ${missing_apt[*]}"
  echo ''
  if $AUTO_YES; then
    :
  else
    read -r -p "Install these now? (will require sudo) [y/N] " response
    if [[ ! "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
      exit
    fi
  fi
  set -x
  sudo apt-get update &&
  sudo apt-get install --yes "${missing_apt[@]}"
}

################################################################################
# Mac (Homebrew)
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

check_mac_deps() {
  missing_brew=()
  check_compiler missing_brew
  # gsl-config is installed alongside GSL and is the standard way to detect it
  if ! $TEST_MODE && command -v gsl-config &>/dev/null; then
    echo "  Found gsl: $(gsl-config --version)"
  else
    echo "  Missing: gsl"
    missing_brew+=(gsl)
  fi
}

install_mac() {
  check_mac_deps
  if [ ${#missing_brew[@]} -eq 0 ]; then
    echo ''
    echo 'All dependencies already installed.'
    exit 0
  fi
  echo ''
  echo "Need to install: ${missing_brew[*]}"
  echo "  Command: brew install ${missing_brew[*]}"
  echo ''
  if $AUTO_YES; then
    :
  else
    read -r -p "Install these now? [y/N] " response
    if [[ ! "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
      exit
    fi
  fi
  install_homebrew
  set -x
  brew install "${missing_brew[@]}"
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
