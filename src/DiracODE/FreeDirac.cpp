#include "DiracODE/FreeDirac.hpp"
#include "Angular/Wigner369j.hpp"
#include "DiracODE/ContinuumState.hpp"
#include "Maths/Grid.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cassert>
#include <cmath>
#include <cstdint>
#include <gsl/gsl_errno.h>
#include <memory>
#include <vector>

namespace DiracODE {

//==============================================================================
FreeDiracParameters::FreeDiracParameters(double en, const Grid &grid,
                                         double alpha) {
  assert(en > 0.0 && "free Dirac waves require a positive (continuum) energy");
  // f -> D sin(kr - l pi/2), while r j_l(kr) -> sin(kr - l pi/2)/k
  k = std::sqrt(en * (2.0 + alpha * alpha * en));
  norm = k * analytic_f_amplitude(en, alpha);
  // c k / (E + m c^2) = k alpha / (2 + alpha^2 en)
  g_ratio = k * alpha / (2.0 + alpha * alpha * en);
  // Store the solution only where the grid resolves the oscillations (at
  // least N_ppw points per wavelength), as solveContinuum does, so that
  // free and distorted waves at the same energy are truncated alike
  const int N_ppw = 10;
  const double dr_max = 2.0 * M_PI / (k * N_ppw);
  max_pt = grid.num_points();
  while (max_pt > 0 && grid.drdu(max_pt - 1) * grid.du() > dr_max) {
    --max_pt;
  }
}

//==============================================================================
std::vector<DiracSpinor> freeDirac(double en, int min_l, int max_l,
                                   std::shared_ptr<const Grid> grid,
                                   double alpha) {
  assert(grid != nullptr);
  assert(min_l >= 0 && max_l >= min_l);

  // Underflow of j_l(kr) at small r and large l is expected (the value is
  // zero to double precision), not an error: keep GSL from aborting
  [[maybe_unused]] static const auto hndl = gsl_set_error_handler_off();

  const FreeDiracParameters fw(en, *grid, alpha);

  // Every kappa with l in [min_l, max_l], in kappa-index order (-1, 1, -2, ..)
  std::vector<DiracSpinor> waves;
  for (std::uint64_t ki = 0;; ++ki) {
    const auto kappa = Angular::kindex_to_kappa(ki);
    const auto l = Angular::l_k(kappa);
    if (l > max_l)
      break;
    if (l < min_l)
      continue;
    waves.emplace_back(0, kappa, grid);
    waves.back().en() = en;
    waves.back().max_pt() = fw.max_pt;
  }

  // j_l(kr) for l = 0..max_l+1 at each point (g needs l+1 when kappa < 0)
  const int l_top = max_l + 1;
  std::vector<double> jl(std::size_t(l_top) + 1);
  for (std::size_t i = 0; i < fw.max_pt; ++i) {
    const auto r = grid->r(i);
    const auto x = fw.k * r;
    for (int l = 0; l <= l_top; ++l) {
      jl[std::size_t(l)] = SphericalBessel::jL(l, x);
    }
    for (auto &Fk : waves) {
      const auto l = std::size_t(Fk.l());
      const auto l_tilde = std::size_t(Angular::l_tilde_k(Fk.kappa()));
      const double sign = Fk.kappa() > 0 ? 1.0 : -1.0;
      Fk.f(i) = fw.norm * r * jl[l];
      Fk.g(i) = sign * fw.norm * fw.g_ratio * r * jl[l_tilde];
    }
  }
  return waves;
}

//==============================================================================
DiracSpinor freeDirac(double en, int kappa, std::shared_ptr<const Grid> grid,
                      double alpha) {
  assert(grid != nullptr);
  [[maybe_unused]] static const auto hndl = gsl_set_error_handler_off();

  const FreeDiracParameters fw(en, *grid, alpha);
  const auto l = Angular::l_k(kappa);
  const auto l_tilde = Angular::l_tilde_k(kappa);
  const double sign = kappa > 0 ? 1.0 : -1.0;

  DiracSpinor Fk(0, kappa, grid);
  Fk.en() = en;
  Fk.max_pt() = fw.max_pt;
  for (std::size_t i = 0; i < fw.max_pt; ++i) {
    const auto r = grid->r(i);
    const auto x = fw.k * r;
    Fk.f(i) = fw.norm * r * SphericalBessel::jL(l, x);
    Fk.g(i) = sign * fw.norm * fw.g_ratio * r * SphericalBessel::jL(l_tilde, x);
  }
  return Fk;
}

} // namespace DiracODE
