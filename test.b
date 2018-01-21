#!/bin/bash
# -std=c++11
g++ -O -Wall -Wpedantic -o test.x test.cpp test.h MAT_matrixAlgebraGSL.h adamsSolveLocalBS.h physicalConstants.h atomInfo.h adamsSolveLocalBS.h -lgsl -lgslcblas -lm &&
./test.x
