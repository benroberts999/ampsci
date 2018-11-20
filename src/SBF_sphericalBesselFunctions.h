#ifndef _SBF_H
#define _SBF_H
#include <gsl/gsl_sf_bessel.h>
#include <cmath>

namespace SBF{

double JL(int L, double x);
double exactGSL_JL(int L, double x);

}

#endif
