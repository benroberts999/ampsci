#!/bin/bash
make clean && make CXX='g++-9 --coverage -g' Build=dev tests &&
./tests [unit] &&
lcov --capture --directory . --output-file coverage.info && 
# don't include test framework itsef in report
lcov --remove coverage.info '*/catch2/*' --output-file coverage.info && 
lcov --remove coverage.info '/*.tests*' --output-file coverage.info && 
# don't include system headers
lcov --remove coverage.info '/usr/*' --output-file coverage.info && 
# don't include Modules - these just run other parts of the code
lcov --remove coverage.info '*/Modules/*' --output-file coverage.info && 
lcov --remove coverage.info '*/DMionisation/*' --output-file coverage.info &&
# TEMPORARY: exclude ExternalField and MBPT - these are well tested,
# but the tests are too slow to run in CI
lcov --remove coverage.info '*/ExternalField/*' --output-file coverage.info && 
lcov --remove coverage.info '*/MBPT/*' --output-file coverage.info && 
lcov --list coverage.info |tee -a cov-info.txt
