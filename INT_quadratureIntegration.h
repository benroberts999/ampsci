#ifndef _INTQUAD_H
#define _INTQUAD_H
#include <vector>

double INT_integrate(std::vector<double> f, std::vector<double> w, double h, int l,
  int m, int nquad=14);

#endif
