#pragma once
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {

void boundState(DiracSpinor &psi, const double en0,
                const std::vector<double> &v, const std::vector<double> &H_mag,
                const double alpha, int log_eps = 14);

void regularAtOrigin(DiracSpinor &phi, const double en,
                     const std::vector<double> &v,
                     const std::vector<double> &H_mag, const double alpha);

void regularAtInfinity(DiracSpinor &phi, const double en,
                       const std::vector<double> &v,
                       const std::vector<double> &H_mag, const double alpha);

void solveContinuum(DiracSpinor &phi, const double en,
                    const std::vector<double> &v, const Grid &ext_grid,
                    const double r_asym0, const double alpha);

//******************************************************************************
DiracSpinor solve_inhomog(const int kappa, const double en,
                          const std::vector<double> &v,
                          const std::vector<double> &H_mag, const double alpha,
                          const DiracSpinor &source);

void solve_inhomog(DiracSpinor &phi, const double en,
                   const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source);

void solve_inhomog(DiracSpinor &phi, DiracSpinor &phi0, DiracSpinor &phiI,
                   const double en, const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source);

} // namespace DiracODE
