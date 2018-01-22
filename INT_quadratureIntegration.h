#ifndef _INTQUAD_H
#define _INTQUAD_H
#include <vector>

double INT_integrate(std::vector<double> f, std::vector<double> w, double h,
  int l=0, int m=0, int nquad=14);

#endif
