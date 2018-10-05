//https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=3j#coupling-coefficients
#ifndef _WIGNER_H
#define _WIGNER_H
#include <gsl/gsl_sf_coupling.h>
#include <math.h>
namespace WIG{

double threej(double j1, double j2, double j3, double m1, double m2, double m3);
double threej_1(int j1, int j2, int j3, int m1, int m2, int m3);
double threej_2(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2,
  int two_m3);

double cg(double j1, double m1, double j2, double m2, double J, double M);
double cg_1(int j1, int m1, int j2, int m2, int J, int M);
double cg_2(int two_j1, int two_m1, int two_j2, int two_m2, int two_J,
  int two_M);

double sixj(double j1, double j2, double j3, double j4, double j5, double j6);
double sixj_1(int j1, int j2, int j3, int j4, int j5, int j6);
double sixj_2(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5,
  int two_j6);

double ninej(double j1, double j2, double j3, double j4, double j5, double j6,
  double j7, double j8, double j9);
double ninej_1(int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8,
  int j9);
double ninej_2(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5,
  int two_j6, int two_j7, int two_j8, int two_j9);

int parity(int la, int lb, int k);

}
#endif
