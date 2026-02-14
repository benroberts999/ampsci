#pragma once
#include "AdamsMoulton.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <complex>
#include <utility>
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {

//! For given complex energy en, solves DE with correct boundary conditions at the origin
void regularAtOrigin_C(DiracSpinor &FaR, DiracSpinor &FaI,
                       const std::complex<double> en,
                       const std::vector<double> &v,
                       const std::vector<double> &H_off_diag,
                       const double alpha);

//! For given complex energy en, solves (local) DE with correct boundary conditions at infinity
void regularAtInfinity_C(DiracSpinor &FaR, DiracSpinor &FaI,
                         const std::complex<double> en,
                         const std::vector<double> &v,
                         const std::vector<double> &H_off_diag,
                         const double alpha);

namespace Internal {
//==============================================================================
//! Matrix which defines Dirac derivative: (dF/dr) = D*F
struct CDiracDerivative
  : AdamsMoulton::DerivativeMatrix<std::size_t, std::complex<double>> {

  CDiracDerivative(const Grid &in_grid, const std::vector<double> &in_v,
                   const int in_k, const std::complex<double> in_en,
                   const double in_alpha,
                   const std::vector<double> &V_off_diag = {});
  const Grid *const pgr;
  const std::vector<double> *const v;
  const std::vector<double> *const Hmag;
  const double zion = 1.0;
  const int k;
  const std::complex<double> en;
  const double alpha, cc;

  std::complex<double> a(std::size_t i) const final;
  std::complex<double> b(std::size_t i) const final;
  std::complex<double> c(std::size_t i) const final;
  std::complex<double> d(std::size_t i) const final;

  CDiracDerivative(const CDiracDerivative &) = delete;
  void operator=(const CDiracDerivative &) = delete;
};

// Solves Dirac equation by integrating outwards from zero.
// Integrates only to 'final' (not inclusive). If final=0, goes to f.size()
// Solution has correct boundary condition at r=0, but not at large r.
void solve_Dirac_outwards_C(std::vector<std::complex<double>> &f,
                            std::vector<std::complex<double>> &g,
                            const CDiracDerivative &Hd, std::size_t final = 0);

// Solves Dirac equation by integrating inwards from 'pinf' to 'ctp'
// Solution has correct boundary condition at large r, but not at small r.
void solve_Dirac_inwards_C(std::vector<std::complex<double>> &f,
                           std::vector<std::complex<double>> &g,
                           const CDiracDerivative &Hd, std::size_t ctp,
                           std::size_t pinf);

} // namespace Internal
} // namespace DiracODE