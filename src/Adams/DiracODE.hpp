#pragma once
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {

void boundState(DiracSpinor &psi, const std::vector<double> &v,
                const Grid &rgrid, const double alpha, int log_dele = 0);

void regularAtOrigin(DiracSpinor &phi, const double en,
                     const std::vector<double> &v, const double alpha);
void regularAtInfinity(DiracSpinor &phi, const double en,
                       const std::vector<double> &v, const double alpha);

int solveContinuum(DiracSpinor &phi, const std::vector<double> &v,
                   const Grid &ext_grid, std::size_t i_asym, double alpha);

//******************************************************************************
void solve_inhomog(DiracSpinor &phi, const double en,
                   const std::vector<double> &v, const double alpha,
                   const DiracSpinor &source);

DiracSpinor solve_inhomog(const int kappa, const double en,
                          const std::vector<double> &v, const double alpha,
                          const DiracSpinor &source);

void solve_inhomog(DiracSpinor &phi, DiracSpinor &phi0, DiracSpinor &phiI,
                   const double en, const std::vector<double> &v,
                   const double alpha, const DiracSpinor &source);

} // namespace DiracODE
