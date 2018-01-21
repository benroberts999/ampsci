#ifndef _SOLVEBS_H
#define _SOLVEBS_H
#include <iostream>
//#include <fstream>
#include <cmath>
//#include <stdio.h>
//#include <time.h>
#include <vector>
#include "MAT_matrixAlgebraGSL.h"
#include <gsl/gsl_linalg.h> //??? needed??
#include "INT_quadratureIntegration.h"


  const int AMO=7; //XXX


#endif

int solveDBS(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int Z, int n, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int &pinf, int &its, double &eps, double alpha);

int outwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int Z, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int ctp, double alpha);

int inwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int ctp, int pinf, double alpha);

int adamsMoulton(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int ngp,
    int ni, int nf, double alpha);

int getAdamsCoefs(std::vector<double> &mia, double &mid, double &miaa);

int getOutwardCoefs(std::vector< std::vector<double> > &oie,
    std::vector<double> &oia, double &oid);
