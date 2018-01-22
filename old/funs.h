#ifndef _FUNS_H
#define _FUNS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include "params.h"
#include <vector>
#include <gsl/gsl_linalg.h>
#endif

void debug(int lineno);
double drdt(int i);
double r(int i);
double dror(int i);
double coul(int i, int Z);
double sphere(int i, int Z, int A, double cR0);
double green(int i, int Z);
double tietz(int i, int Z);
double TFapprox(int i, int Z);
double integrate(double *f, int l, int m);
int diff(double *f, double *deriv);
double diracen(int z, int n, int k);


extern "C" void dgetrf_(int*, int*, double*, int*, int*, int*);
extern "C" void dgetri_(int*, double*, int*, int*, double*, int*, int*);
//int invertmat(double (*matrix)[amo2], double (*inverse)[amo2], int dim);
int invertmat(double matrix[amo2][amo2], double inverse[amo2][amo2], int dim);

int invertMatrix(double inmat[amo2][amo2], double outmat[amo2][amo2]);
