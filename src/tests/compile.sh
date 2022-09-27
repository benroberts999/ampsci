#!/bin/bash
set -e

# openMP doesn't seem to work with these:
compilers_list_no_omp=("clang++-6.0" "clang++-7" "clang++-8" "clang++-9" "clang++-11")
# openMP is supported with these compilers
compilers_list=("clang++-10" "g++-7" "g++-8" "g++-9" "g++-10" "g++-11")

for cxx in ${compilers_list_no_omp[@]}; do
  if [ -x "$(command -v "$cxx")" ]; then
    echo $cxx
    make clean
    make Build=dev CXX=$cxx' -Werror' UseOpenMP=no
  fi
done

for cxx in ${compilers_list[@]}; do
  if [ -x "$(command -v "$cxx")" ]; then
    echo $cxx
    make clean
    make Build=dev CXX=$cxx' -Werror' UseOpenMP=yes
  fi
done

./tests [~Slow] |tee -a examples.out
./ampsci Xe |tee -a examples.out
./ampsci ./doc/examples/ampsci.in |tee -a examples.out
./ampsci ./doc/examples/Cs_testBasis.in |tee -a examples.out
./ampsci ./doc/examples/Cs_basicExample.in |tee -a examples.out
./ampsci ./doc/examples/Cs_correlationsExample.in |tee -a examples.out
echo "Success!"