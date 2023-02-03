#!/bin/bash
make clean && make CXX='g++-11' CARGS='--coverage -g' LARGS='--coverage -g' Build=dev tests &&
./tests [unit] &&
lcov --capture --directory . --output-file coverage.info &&
# # don't include test framework itsef in report
# # don't include system headers
# # don't include Modules - these just run other parts of the code
# # TEMPORARY: exclude some from MBPT - these are well tested,
# # but the tests are too slow to run in CI
lcov --remove coverage.info \
  '*/catch2/*' \
  '*/fmt/*' \
  '/*.tests*' \
  '*/version/*' \
  '/usr/*' \
  '*/Modules/*' \
  '*/DMionisation/*' \
  '*/Physics/periodicTable.*' \
  '*/MBPT/FeynmanSigma.*' \
  '*/MBPT/CorrelationPotential.*' \
  '*/MBPT/Ladder.*' \
  '*/qip/Widgets.*' \
--output-file coverage.info &&
lcov --list coverage.info |tee -a cov-info.txt
