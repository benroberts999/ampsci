#ifndef _MATRIXALG_H
#define _MATRIXALG_H
#include <gsl/gsl_linalg.h>
#include <vector>
#include <cmath>  //?? needed?

namespace MAT{

int invertMatrix(std::vector< std::vector<double> > inmat,
                        std::vector< std::vector<double> > &outmat);
int invertMatrix(std::vector< std::vector<double> > &inmat);
int invertMatrix(std::vector< std::vector<float> > inmat,
                        std::vector< std::vector<float> > &outmat);
int invertMatrix(std::vector< std::vector<float> > &inmat);

double calcDeterminant(std::vector< std::vector<double> > inmat);
double calcDeterminant(std::vector< std::vector<float> > inmat);

int linsolve(std::vector< std::vector<double> > inmat,
                 std::vector<double> invec, std::vector<double> &outvec);

}

#endif
