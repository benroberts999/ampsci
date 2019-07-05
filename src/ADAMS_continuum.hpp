#pragma once
#include <vector>
class DiracSpinor;
class Grid;

namespace ADAMS {

double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
                    double y3);

int solveContinuum(DiracSpinor &phi, const std::vector<double> &v,
                   const Grid &ext_grid, std::size_t i_asym, double alpha);

double findSineAmplitude(std::vector<double> &pc, const std::vector<double> &rc,
                         std::size_t NGPc, std::size_t i_asym);

std::size_t findAsymptoticRegion(std::vector<double> &pc,
                                 const std::vector<double> &rc,
                                 std::size_t NGPb, std::size_t NGPc,
                                 std::size_t i_asym);

} // namespace ADAMS
