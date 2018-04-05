#ifndef _SOLVECNTM_H
#define _SOLVECNTM_H
#include "adamsSolveLocalBS.h"


double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
  double y3);

int solveContinuum(std::vector<double> &p, std::vector<double> &q, double en,
    std::vector<double> v, int Z, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    double alpha);


#endif
