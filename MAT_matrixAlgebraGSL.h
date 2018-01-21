#ifndef _MATRIXALG_H
#define _MATRIXALG_H
#include <gsl/gsl_linalg.h>
#include <vector>
#include <cmath>  //?? needed?


int MAT_invertMatrix(std::vector< std::vector<double> > inmat,
                        std::vector< std::vector<double> > &outmat);
int MAT_invertMatrix(std::vector< std::vector<double> > &inmat);
int MAT_invertMatrix(std::vector< std::vector<float> > inmat,
                        std::vector< std::vector<float> > &outmat);
int MAT_invertMatrix(std::vector< std::vector<float> > &inmat);

double MAT_calcDeterminant(std::vector< std::vector<double> > inmat);
double MAT_calcDeterminant(std::vector< std::vector<float> > inmat);

int MAT_linsolve(std::vector< std::vector<double> > inmat,
                 std::vector<double> invec, std::vector<double> &outvec);

#endif
