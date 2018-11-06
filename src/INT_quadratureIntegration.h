#ifndef _INTQUAD_H
#define _INTQUAD_H
#include <vector>
#include <iostream>

namespace INT{

double integrate(std::vector<double> f, std::vector<double> w, double h,
  int l=0, int m=0, int nquad=14);

int diff(std::vector<double> f, std::vector<double> drdt, double h,
      std::vector<double> &deriv);

double integrate4(
        std::vector<double> &f1,
        std::vector<double> &f2,
        std::vector<double> &f3,
        std::vector<double> &f4,
        double h=1, int l=0, int m=0, int a_start=1, int a_end=1);
double integrate3(
        std::vector<double> &f1,
        std::vector<double> &f2,
        std::vector<double> &f3,
        double h=1, int l=0, int m=0, int a_start=1, int a_end=1);
double integrate2(std::vector<double> &f1,std::vector<double> &f2,
        double h=1, int l=0, int m=0,
      int a_start=1, int a_end=1);
double integrate1(std::vector<double> &f1, double h=1, int l=0, int m=0,
   int a_start=1, int a_end=1);

}

#endif
