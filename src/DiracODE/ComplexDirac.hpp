#pragma once
#include "AdamsMoulton.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <complex>
#include <utility>
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {

//! For given complex energy en, solves Dirac equation with correct boundary conditions at the origin.
void regularAtOrigin_C(DiracSpinor &FaR, DiracSpinor &FaI,
                       const std::complex<double> en,
                       const std::vector<double> &v,
                       const std::vector<double> &H_off_diag,
                       const double alpha);

//! For given complex energy en, solves Dirac equation with correct boundary conditions at infinity.
void regularAtInfinity_C(DiracSpinor &FaR, DiracSpinor &FaI,
                         const std::complex<double> en,
                         const std::vector<double> &v,
                         const std::vector<double> &H_off_diag,
                         const double alpha);

namespace Internal {
//==============================================================================
/*!
  @brief Complex-energy Dirac derivative matrix: dF/du = D(u)*F(u).
  @details
  Complex-energy analogue of @ref DiracDerivative. Accepts a complex orbital
  energy en, returning complex matrix elements. See @ref DiracDerivative for
  the form of D; the only difference here is Y = std::complex<double>.

  @note Non-copyable; stores raw pointers to grid and potential arrays.
*/
struct CDiracDerivative
  : AdamsMoulton::DerivativeMatrix<std::size_t, std::complex<double>> {

  /*!
    @brief Constructs the complex Dirac derivative matrix.
    @param in_grid     Radial grid.
    @param in_v        Local potential v(r).
    @param in_k        Orbital kappa quantum number.
    @param in_en       Complex orbital energy.
    @param in_alpha    Fine-structure constant.
    @param V_off_diag  Optional off-diagonal (magnetic) potential. If empty, treated as zero.
  */
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

  //! D matrix elements (see @ref DiracDerivative for definitions); index i is grid point.
  std::complex<double> a(std::size_t i) const final;
  std::complex<double> b(std::size_t i) const final;
  std::complex<double> c(std::size_t i) const final;
  std::complex<double> d(std::size_t i) const final;

  CDiracDerivative(const CDiracDerivative &) = delete;
  void operator=(const CDiracDerivative &) = delete;
};

/*!
  @brief Integrates the complex Dirac equation outwards from the origin.
  @details
  Integrates up to index @p final (not inclusive); if final=0, integrates to
  f.size(). Solution satisfies the boundary condition at r=0 but not at large r.
*/
void solve_Dirac_outwards_C(std::vector<std::complex<double>> &f,
                            std::vector<std::complex<double>> &g,
                            const CDiracDerivative &Hd, std::size_t final = 0);

//! Integrates the complex Dirac equation inwards from pinf to ctp.
//! Solution satisfies the boundary condition at large r but not at the origin.
void solve_Dirac_inwards_C(std::vector<std::complex<double>> &f,
                           std::vector<std::complex<double>> &g,
                           const CDiracDerivative &Hd, std::size_t ctp,
                           std::size_t pinf);

} // namespace Internal
} // namespace DiracODE