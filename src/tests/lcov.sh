#!/bin/bash
make clean && make CXX='g++-9 --coverage -g' Build=dev tests &&
./tests [unit] &&
lcov --capture --directory . --output-file coverage.info &&
# # don't include test framework itsef in report
# # don't include system headers
# # don't include Modules - these just run other parts of the code
# # TEMPORARY: exclude ExternalField and MBPT - these are well tested,
# # but the tests are too slow to run in CI
lcov --remove coverage.info \
  '*/catch2/*' \
  '/*.tests*' \
  '*/version/*' \
  '/usr/*' \
  '*/Modules/*' \
  '*/DMionisation/*' \
  '*/ExternalField/*' \
  '*/MBPT/*' \
--output-file coverage.info &&
lcov --list coverage.info |tee -a cov-info.txt
