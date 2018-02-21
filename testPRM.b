#!/bin/bash
g++ -std=c++11 -O -Wall -Wpedantic -fopenmp -o testPRM.x testPRM.cpp ElectronOrbitals.cpp \
MAT_matrixAlgebraGSL.cpp adamsSolveLocalBS.cpp \
INT_quadratureIntegration.cpp PRM_parametricPotentials.cpp -lgsl -lgslcblas -lm &&
./testPRM.x
