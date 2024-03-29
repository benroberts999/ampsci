#!/bin/bash
set -e

# openMP doesn't seem to work with these:
compilers_list=("clang++-6.0" "clang++-7" "clang++-8" "clang++-9" "clang++-10" "clang++-11" "clang++-12" "clang++-13" "clang++-14" "clang++-15" "g++-7" "g++-8" "g++-9" "g++-10" "g++-11" "g++-12 g++-13")
# openMP is supported with these compilers
compilers_list_omp=("clang++-12" "clang++-13" "clang++-14" "g++-7" "g++-8" "g++-9" "g++-10" "g++-11" "g++-12 g++-13")

for cxx in ${compilers_list[@]}; do
  if [ -x "$(command -v "$cxx")" ]; then
    echo $cxx
    # make clean
    make Build=dev CXX=$cxx' -Werror' UseOpenMP=no
  fi
done

for cxx in ${compilers_list_omp[@]}; do
  if [ -x "$(command -v "$cxx")" ]; then
    echo $cxx
    # make clean
    make Build=dev CXX=$cxx' -Werror' UseOpenMP=yes
  fi
done

./tests [unit] |tee -a examples.out
./ampsci Xe |tee -a examples.out
./ampsci ./doc/examples/ampsci.in |tee -a examples.out
./ampsci ./doc/examples/Cs_testBasis.in |tee -a examples.out
./ampsci ./doc/examples/Cs_basicExample.in |tee -a examples.out
./ampsci ./doc/examples/Cs_correlationsExample.in |tee -a examples.out
echo "Success!"
