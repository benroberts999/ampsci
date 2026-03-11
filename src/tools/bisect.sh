#!/usr/bin/env bash

git bisect start
git bisect bad  86233288
git bisect good 4e250f2e

# Account for change in makefile format
git bisect run bash -lc '
  rm Makefile
  cp ./doc/Makefile .              2>/dev/null || :
  cp ./doc/examples/Makefile .     2>/dev/null || :
  make clean
  make tests -j12 MODE=dev Build=dev
  ./tests "MBPT: Correlation Potential: Sigma2"
'