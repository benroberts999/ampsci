#!/bin/bash
for file in $(find src -name '*.cpp' -or -name '*.hpp'); do
  while true; do
    rm -f obj/junk_abcde123_gnu.o
    rm -f junk_abcde123_clang.o
    echo ""
    echo "Compiling file: "$file
    echo ""
    echo "with GNU..."
    g++ -std=c++14 -c -Wno-unknown-pragmas -Wpedantic -Wall -Wextra -Wdouble-promotion -Wconversion -Wlogical-op $file -o obj/junk_abcde123_gnu.o -I./src/
    if test -f "junk_abcde123_gnu.o"; then
      echo "done"
    fi
    rm -f obj/junk_abcde123_gnu.o
  #
    echo ""
    echo "with CLANG..."
    clang++ -std=c++14 -c -Wno-unknown-pragmas -Wpedantic -Wall -Wextra -Wdouble-promotion -Wconversion -Wno-sign-conversion -Wheader-hygiene $file -o obj/junk_abcde123_clang.o -I./src/
    if test -f "junk_abcde123_clang.o"; then
      echo "done clang"
    fi
    rm -f obj/junk_abcde123_clang.o
  #
    echo ""
    echo "clang tidy:"
    clang-tidy $file -checks=clang-analyzer-*,openmp-*,performance-* -- -I./src/
    echo "done"
  #
    echo ""
    echo "Finished: "$file
    echo "Press enter to continue, or any other key to re-check $file:"
    read -p "... " redo
    if [ "$redo" == "" ]; then
      break;
    fi
  done
done