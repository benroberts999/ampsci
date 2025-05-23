#!/bin/bash
make clean && make CXX='g++-11' CARGS='--coverage -g' LARGS='--coverage -g' Build=dev tests &&
rm -f ./coverage.info &&
lcov --gcov-tool gcov-11 --zerocounters --directory . &&
./tests [unit] &&
lcov --gcov-tool gcov-11 --capture --directory . --output-file coverage.info &&
# # don't include test framework itsef in report
# # don't include system headers
# # don't include Modules - these just run other parts of the code
lcov --gcov-tool gcov-11 --remove coverage.info \
  '*/catch2/*' \
  '*/fmt/*' \
  '/*.tests*' \
  '*/version/*' \
  '/usr/*' \
  '*/Modules/*' \
  '*/Kionisation/*' \
  '*/Physics/periodicTable.*' \
  '*/DiracODE/ComplexDirac.*' \
  '*/MBPT/Ladder.*' \
  '*/qip/Widgets.*' \
--output-file coverage.info &&
lcov --list coverage.info |tee -a cov-info.txt

# https://lcov-viewer.netlify.app/
# https://github.com/eugenezinovyev/lcov-viewer