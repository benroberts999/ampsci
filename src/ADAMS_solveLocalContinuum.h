#pragma once
#include <vector>
namespace ADAMS {

double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
                    double y3);

int solveContinuum(std::vector<double> &f, std::vector<double> &g, double en,
                   const std::vector<double> &v, int ka,
                   const std::vector<double> &rc,
                   const std::vector<double> &drdt, double h, int NGPb,
                   int NGPc, int i_asym, double alpha);

double findSineAmplitude(std::vector<double> &pc, const std::vector<double> &rc,
                         int NGPc, int i_asym);

int findAsymptoticRegion(std::vector<double> &pc, const std::vector<double> &rc,
                         int NGPb, int NGPc, int i_asym);

} // namespace ADAMS
