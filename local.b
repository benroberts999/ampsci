#!/bin/bash
g++ -O -Wall -Wpedantic -o local.x local.cpp solvebs.cpp solvebs.h funs.cpp funs.h params.h -lgsl -lgslcblas -lm -llapack -lblas &&
./local.x
