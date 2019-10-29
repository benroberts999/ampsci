#pragma once
#include "DiracODE.hpp"
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {
namespace Adams {

double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
                    double y3);

double findSineAmplitude(std::vector<double> &pc, const std::vector<double> &rc,
                         std::size_t num_pointsc, std::size_t i_asym);

std::size_t findAsymptoticRegion(std::vector<double> &pc,
                                 const std::vector<double> &rc,
                                 std::size_t num_pointsb,
                                 std::size_t num_pointsc, std::size_t i_asym);

} // namespace Adams
} // namespace DiracODE
