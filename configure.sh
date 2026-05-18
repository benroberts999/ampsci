#!/bin/bash

echo ''
echo 'This will create a Makefile to compile ampsci using the default options.'
echo 'You may need to update some options in the Makefile, depending on your system.'
echo 'You may need to run install-dependencies.sh first, to install dependencies.'
echo ''

if [[ "$1" != "--yes" && "$1" != "-y" ]]; then
  read -r -p "This script will overwrite existing Makefile. Do you wish to continue? [y/N] " response
  if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
  then
    echo ""
  else
    exit 0
  fi
fi

# Figure out which operating system is used:
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    *)          machine="UNKNOWN:${unameOut}"
esac
echo "Configuring for : ${machine}"
if [[ "${machine}" == "UNKNOWN:${unameOut}" ]]; then
  echo 'Unsupported system for setup - requires manual setup'
  exit 1
fi

# Copy template Makefile and example input file
cp ./doc/Makefile ./ || { echo "Error: doc/Makefile not found"; exit 1; }
cp -n ./doc/examples/ampsci.in ./ampsci.in

################################################################################
# Compiler detection
# Search for the highest available versioned g++ (e.g. g++-14), then fall back
# to versioned clang++, then unversioned g++/clang++.
# Using versioned names ensures we get a real GCC on Mac, where plain g++ is
# Apple Clang in disguise.
################################################################################

# Minimum supported compiler versions:
min_gcc_ver=7
min_clang_ver=6

found_compiler=0
is_clang=0

# Search for highest available versioned GCC (g++-N), counting down from max:
gccver=30
while (( gccver >= min_gcc_ver )); do
  if [[ -x "$(command -v "g++-${gccver}")" ]]; then
    found_compiler=1
    cxx=$(command -v "g++-${gccver}")
    echo "Compiler Found  : ${cxx}"
    sed -i.bak "s@^CXX[[:space:]]*=.*@CXX = ${cxx}@g" Makefile
    break
  fi
  (( gccver-- ))
done

# If no versioned GCC, search for highest available versioned clang++:
if (( found_compiler == 0 )); then
  clangver=30
  while (( clangver >= min_clang_ver )); do
    if [[ -x "$(command -v "clang++-${clangver}")" ]]; then
      found_compiler=1
      is_clang=1
      cxx=$(command -v "clang++-${clangver}")
      echo "Compiler Found  : ${cxx}"
      sed -i.bak "s@^CXX[[:space:]]*=.*@CXX = ${cxx}@g" Makefile
      break
    fi
    (( clangver-- ))
  done
fi

# Last resort: try unversioned g++ or clang++
if (( found_compiler == 0 )); then
  if [[ -x "$(command -v "g++")" ]]; then
    found_compiler=1
    cxx="g++"
    echo "Compiler Found  : ${cxx}"
    sed -i.bak "s@^CXX[[:space:]]*=.*@CXX = ${cxx}@g" Makefile
  elif [[ -x "$(command -v "clang++")" ]]; then
    found_compiler=1
    is_clang=1
    cxx="clang++"
    echo "Compiler Found  : ${cxx}"
    sed -i.bak "s@^CXX[[:space:]]*=.*@CXX = ${cxx}@g" Makefile
  fi
fi

if (( found_compiler == 0 )); then
  echo "Compatible gcc or clang compiler not found."
  echo "Run: install-dependencies.sh, or"
  if [[ "${machine}" == 'Linux' ]]; then
    echo "Run: sudo apt install g++"
  elif [[ "${machine}" == 'Mac' ]]; then
    echo "Run: brew install gcc"
  fi
  echo "Or proceed with manual compilation for other compilers"
  exit 1
fi

# Check compiler version meets minimum requirement.
# -dumpversion prints the major version number on modern compilers.
cxx_ver=$("$cxx" -dumpversion 2>/dev/null | cut -d. -f1)
if (( is_clang == 1 )); then
  min_ver=$min_clang_ver
else
  min_ver=$min_gcc_ver
fi
if [[ -n "$cxx_ver" ]] && (( cxx_ver < min_ver )); then
  echo "Compiler version: ${cxx_ver} -- WARNING: below minimum (${min_ver}), may not compile"
else
  echo "Compiler version: ${cxx_ver} (ok)"
fi

################################################################################
# GSL library detection
# gsl-config is the standard tool installed alongside GSL; it knows the prefix
# regardless of whether GSL is in a system path or a Homebrew/custom location.
################################################################################
if [[ -x "$(command -v "gsl-config")" ]]; then
  gslver=$(gsl-config --version)
  gslpath=$(gsl-config --prefix)
  gsllibs=$(gsl-config --libs)
  echo "GSL version     : ${gslver}"
  echo "GSL path        : ${gslpath}"
  echo "GSL flags       : ${gsllibs}"
  sed -i.bak "s@^GSL_PATH.*@GSL_PATH ?= ${gslpath}@g" Makefile
else
  echo "GSL libraries not installed"
  echo "Run: install-dependencies.sh, or"
  if [[ "${machine}" == 'Linux' ]]; then
    echo "Run: sudo apt install libgsl-dev"
  elif [[ "${machine}" == 'Mac' ]]; then
    echo "Run: brew install gsl"
  fi
  exit 1
fi

################################################################################
# LAPACK/BLAS linker flags
# On Mac: use Apple's Accelerate framework (built-in, no install needed).
# On Linux: try pkg-config first, then detect OpenBLAS (used by foss and other
# HPC toolchains), then fall back to standard -llapack -lblas.
# EBROOTOPENBLAS is set by EasyBuild when OpenBLAS module is loaded.
# GSL flags come from gsl-config (already run above).
################################################################################
if [[ "${machine}" == "Mac" ]]; then
  ldlibs="${gsllibs} -framework Accelerate"
  echo "LDLIBS          : ${ldlibs}"
  sed -i.bak "s@^LDLIBS ?=.*@LDLIBS ?= ${ldlibs}@" Makefile
else
  blas_libs=""
  if command -v pkg-config &>/dev/null; then
    blas_libs=$(pkg-config --libs lapack blas 2>/dev/null)
  fi

  if [[ -z "$blas_libs" ]]; then
    # pkg-config failed -- detect manually
    # OpenBLAS includes LAPACK; preferred over separate -llapack -lblas
    if [[ -n "$EBROOTOPENBLAS" ]] || \
       echo "int main(){return 0;}" | "$cxx" -lopenblas -x c++ - -o /dev/null 2>/dev/null; then
      blas_libs="-lopenblas"
      # foss/EasyBuild OpenBLAS builds often need -lgfortran for Fortran runtime
      if [[ -n "$EBROOTGCC" || -n "$EBROOTGCCCORE" ]]; then
        blas_libs="${blas_libs} -lgfortran"
      fi
    else
      blas_libs="-llapack -lblas"
    fi
  fi

  ldlibs="${gsllibs} ${blas_libs}"
  echo "LDLIBS          : ${ldlibs}"
  sed -i.bak "s@^LDLIBS ?=.*@LDLIBS ?= ${ldlibs}@" Makefile
fi

################################################################################
# Parallel build jobs
# nproc (Linux) and sysctl hw.logicalcpu (Mac) both report logical CPU count.
################################################################################
if command -v nproc &>/dev/null; then
  ncpus=$(nproc)
elif command -v sysctl &>/dev/null; then
  ncpus=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)
else
  ncpus=1
fi
(( ncpus > 8 )) && ncpus=8
echo "Parallel jobs   : ${ncpus}"
sed -i.bak "s@^DEFAULT_MFLAGS.*@DEFAULT_MFLAGS = -j${ncpus}@" Makefile

################################################################################
# Build mode: release on main branch, dev on all other branches.
# Dev mode enables extra warnings and disables some optimisations, which is
# more useful when actively editing code.
################################################################################
if command -v git &>/dev/null && git rev-parse --is-inside-work-tree &>/dev/null; then
  branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null)
  if [[ "$branch" == "main" ]]; then
    build_mode="release"
  else
    build_mode="dev"
  fi
  echo "Git branch      : ${branch}"
  echo "Build mode      : ${build_mode}"
  sed -i.bak "s@^MODE ?=.*@MODE ?= ${build_mode}@" Makefile
fi

################################################################################
# OpenMP detection
# Try flags in order; the first one that compiles a trivial program is used.
# -fopenmp=libomp is needed for clang on Mac with Homebrew libomp.
# If none work, OMPLIB is left blank (OpenMP disabled).
################################################################################
ompflag=""
for flag in "-fopenmp" "-fopenmp=libomp" "-fopenmp=libgomp"; do
  if echo "int main(){return 0;}" | "$cxx" $flag -x c++ - -o /dev/null >/dev/null 2>&1; then
    ompflag="$flag"
    break
  fi
done
if [[ -n "$ompflag" ]]; then
  echo "OpenMP          : ${ompflag}"
  sed -i.bak "s@^OMPLIB.*@OMPLIB ?= ${ompflag}@" Makefile
else
  echo "OpenMP          : not supported"
  sed -i.bak "s@^OMPLIB.*@OMPLIB ?=@" Makefile
fi

################################################################################
# Done -- clean up sed backup and report
################################################################################

echo ''
rm -fv Makefile.bak

echo ''
echo "Makefile has been configured"
echo "Attempt to compile ampsci using Make:"
echo "$ make"
echo ''
