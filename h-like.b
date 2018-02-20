#!/bin/bash
g++ -std=c++11 -O -Wall -Wpedantic -o h-like.x h-like.cpp ElectronOrbitals.cpp \
MAT_matrixAlgebraGSL.cpp adamsSolveLocalBS.cpp \
INT_quadratureIntegration.cpp -lgsl -lgslcblas -lm &&
./h-like.x
