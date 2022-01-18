#!/bin/bash
set -e
make clean && make CXX=clang++-6.0 UseOpenMP=no Build=release
make clean && make CXX=clang++-10 UseOpenMP=yes Build=dev
./unitTests quick
make clean && make CXX=g++-7 UseOpenMP=yes Build=release
make clean && make CXX='g++-10 -Werror' UseOpenMP=yes Build=dev
./unitTests quick
./ampsci Xe
./ampsci Br
./ampsci ./doc/examples/ampsci.in
./ampsci ./doc/examples/Cs_testBasis.in
./ampsci ./doc/examples/Cs_basicExample.in
./ampsci ./doc/examples/Cs_correlationsExample.in
echo "Success!"
