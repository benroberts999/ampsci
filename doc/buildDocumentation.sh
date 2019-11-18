#!/bin/bash
file=./doc/03-diracSCAS_code.md
src=./src
echo "# Code documentation:" > $file
echo "" >> $file
echo "(very unfinished..)" >> $file
echo "" >> $file

cat $src/Dirac/Dirac.md >> $file
echo "" >> $file
echo "---" >> $file
echo "" >> $file
cat $src/Adams/Adams.md >> $file
echo "" >> $file
echo "---" >> $file
echo "" >> $file
