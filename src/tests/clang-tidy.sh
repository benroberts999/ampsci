#!/bin/bash
ask_each=false
for file in $(find src -name '*.cpp' -or -name '*.hpp'); do
  while true; do
    echo ""
    echo "clang tidy:"
    PathForGSL=/usr/local/opt/gnu-scientific-library
    clang-tidy $file -extra-arg=-std=c++17 -checks=clang-analyzer-*,openmp-*,performance-*,portability-*,clang-analyzer-* -- -I./src/ -I$PathForGSL/include/
# modernize-*, cppcoreguidelines-*,readability-*,llvm-*,
    echo "done"
  #
    echo ""
    echo "Finished: "$file
    if $ask_each; then
      echo "Press enter to continue, or any other key to re-check $file:"
      read -p "... " redo
      if [ "$redo" == "" ]; then
        break;
      fi
    else
      break;
    fi
  done
done
