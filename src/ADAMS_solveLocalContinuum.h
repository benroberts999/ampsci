#pragma once
#include <vector>
namespace ADAMS {

double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
                    double y3);

int solveContinuum(std::vector<double> &f, std::vector<double> &g, double en,
                   const std::vector<double> &v, int ka,
                   const std::vector<double> &rc,
                   const std::vector<double> &drdt, double h, size_t NGPb,
                   size_t NGPc, size_t i_asym, double alpha);

double findSineAmplitude(std::vector<double> &pc, const std::vector<double> &rc,
                         size_t NGPc, size_t i_asym);

size_t findAsymptoticRegion(std::vector<double> &pc,
                            const std::vector<double> &rc, size_t NGPb,
                            size_t NGPc, size_t i_asym);

} // namespace ADAMS
