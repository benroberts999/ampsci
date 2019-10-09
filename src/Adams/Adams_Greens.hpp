#pragma once
#include <vector>

class DiracSpinor;

namespace Adams {
//******************************************************************************
void solve_inhomog(DiracSpinor &phi, const double en,
                   const std::vector<double> &v, const double alpha,
                   const DiracSpinor &source);

//------------------------------------------------------------------------------
void solve_inhomog(DiracSpinor &phi, DiracSpinor &phi0, DiracSpinor &phiI,
                   const double en, const std::vector<double> &v,
                   const double alpha, const DiracSpinor &source);

//******************************************************************************
void GreenSolution(DiracSpinor &phi, const DiracSpinor &phiI,
                   const DiracSpinor &phi0, const double alpha,
                   const DiracSpinor &Sr);
} // namespace Adams
