#!/bin/bash
g++ -O -Wall -Wpedantic -o test.x test.cpp test.h -lgsl -lgslcblas -lm -llapack -lblas &&
./local.x
