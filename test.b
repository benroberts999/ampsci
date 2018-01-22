#!/bin/bash
# -std=c++11
g++ -std=c++11 -O -Wall -Wpedantic -o test.x test.cpp ElectronOrbitals.cpp ElectronOrbitals.h MAT_matrixAlgebraGSL.cpp MAT_matrixAlgebraGSL.h adamsSolveLocalBS.cpp adamsSolveLocalBS.h physicalConstants.h atomInfo.h INT_quadratureIntegration.cpp INT_quadratureIntegration.h -lgsl -lgslcblas -lm &&
./test.x
