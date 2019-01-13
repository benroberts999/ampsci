#pragma once
#include <vector>

namespace INT{

double integrate(const std::vector<double> &f, const std::vector<double> &w,
  double h, int l=0, int m=0, int nquad=14);

int diff(const std::vector<double> &f, const std::vector<double> &drdt,
      double h, std::vector<double> &deriv, int order=0);

double integrate4(
        const std::vector<double> &f1,
        const std::vector<double> &f2,
        const std::vector<double> &f3,
        const std::vector<double> &f4,
        double h=1, int l=0, int m=0, int a_start=1, int a_end=1);
double integrate3(
        const std::vector<double> &f1,
        const std::vector<double> &f2,
        const std::vector<double> &f3,
        double h=1, int l=0, int m=0, int a_start=1, int a_end=1);
double integrate2(const std::vector<double> &f1,const std::vector<double> &f2,
        double h=1, int l=0, int m=0,
      int a_start=1, int a_end=1);
double integrate1(const std::vector<double> &f1, double h=1, int l=0, int m=0,
   int a_start=1, int a_end=1);

}
