#pragma once
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {

//! @brief For given energy en (en > 0), solves (local) DE for continuum state
//! (with energy normalisation).
/*! @details
ext_grid is an 'extended grid' (see Grid class); needed since we need to solve
to large r to enforce boundary conditions (especially for large en). r_asym0 is
initial guess for 'asymptotic region'. ext_grid must extend past r_asym0.
*/
void solveContinuum(DiracSpinor &Fa, const double en,
                    const std::vector<double> &v, const Grid &ext_grid,
                    const double r_asym0, const double alpha,
                    const DiracSpinor *const VxFa = nullptr,
                    const DiracSpinor *const Fa0 = nullptr);

namespace Adams {

double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
                    double y3);

double findSineAmplitude(const std::vector<double> &pc,
                         const std::vector<double> &rc, std::size_t num_pointsc,
                         std::size_t i_asym);

std::size_t findAsymptoticRegion(const std::vector<double> &pc,
                                 const std::vector<double> &rc,
                                 std::size_t num_pointsb,
                                 std::size_t num_pointsc, std::size_t i_asym);

} // namespace Adams
} // namespace DiracODE
