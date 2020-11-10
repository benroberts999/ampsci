make clean && make CXX='clang++-6.0 -Werror' &&
#make clean && make CXX='clang++-7 -Werror' &&
make clean && make CXX='clang++-9 -Werror' &&
# g++7 still complains unused vars for structured bindings
make clean && make CXX='g++-7 -Wno-unused-variable' &&
#make clean && make CXX='g++-8 -Werror' &&
make clean && make CXX='g++-9 -Werror' &&
./ampsci ./doc/examples/ampsci.in &&
./ampsci ./doc/examples/Cs_basicExample.in &&
./ampsci ./doc/examples/Cs_correlationsExample.in
