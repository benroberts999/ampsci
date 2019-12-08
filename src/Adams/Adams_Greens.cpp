#include "Adams_Greens.hpp"
#include "Adams_bound.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "DiracODE.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include <vector>
/*

Solve inhomogenous Dirac equation:
(H_0 + v - e)phi = S

S (source) is spinor

H_0 includes v_nuc, but not v_el
Usually, V = V_dir, but can be any LOCAL potential.

XXX Update to accept operators?

*/

namespace DiracODE {

DiracSpinor solve_inhomog(const int kappa, const double en,
                          const std::vector<double> &v, const double alpha,
                          const DiracSpinor &source) {
  auto sp = SafeProfiler::profile(__func__, "0");
  auto phi = DiracSpinor(0, kappa, *source.p_rgrid);
  solve_inhomog(phi, en, v, alpha, source);
  return phi;
}

//******************************************************************************
void solve_inhomog(DiracSpinor &phi, const double en,
                   const std::vector<double> &v, const double alpha,
                   const DiracSpinor &source)
// NOTE: returns NON-normalised function!
{
  auto sp = SafeProfiler::profile(__func__, "a");
  auto phi0 = DiracSpinor(phi.n, phi.k, *phi.p_rgrid);
  auto phiI = DiracSpinor(phi.n, phi.k, *phi.p_rgrid);
  solve_inhomog(phi, phi0, phiI, en, v, alpha, source);
}
//------------------------------------------------------------------------------
void solve_inhomog(DiracSpinor &phi, DiracSpinor &phi0, DiracSpinor &phiI,
                   const double en, const std::vector<double> &v,
                   const double alpha, const DiracSpinor &source)
// Overload of the above. Faster, since doesn't need to allocate for phi0 and
// phiI
{
  auto sp = SafeProfiler::profile(__func__, "b");
  regularAtOrigin(phi0, en, v, alpha);
  regularAtInfinity(phiI, en, v, alpha);
  Adams::GreenSolution(phi, phiI, phi0, alpha, source);
}

namespace Adams {
//******************************************************************************
void GreenSolution(DiracSpinor &phi, const DiracSpinor &phiI,
                   const DiracSpinor &phi0, const double alpha,
                   const DiracSpinor &Sr) {
  auto sp = SafeProfiler::profile(__func__);
  // Wronskian. Note: in current method: only need the SIGN
  int pp = int(0.65 * double(phiI.pinf));
  auto w2 = (phiI.f[pp] * phi0.g[pp] - phi0.f[pp] * phiI.g[pp]);
  // XXX test making this non-constant over r?

  // save typing:
  const auto &gr = *phi.p_rgrid;
  constexpr auto ztr = NumCalc::zero_to_r;
  constexpr auto rti = NumCalc::r_to_inf;

  phi *= 0.0;
  NumCalc::additivePIntegral<ztr>(phi.f, phiI.f, phi0.f, Sr.f, gr, phiI.pinf);
  NumCalc::additivePIntegral<ztr>(phi.f, phiI.f, phi0.g, Sr.g, gr, phiI.pinf);
  NumCalc::additivePIntegral<rti>(phi.f, phi0.f, phiI.f, Sr.f, gr, phiI.pinf);
  NumCalc::additivePIntegral<rti>(phi.f, phi0.f, phiI.g, Sr.g, gr, phiI.pinf);
  NumCalc::additivePIntegral<ztr>(phi.g, phiI.g, phi0.f, Sr.f, gr, phiI.pinf);
  NumCalc::additivePIntegral<ztr>(phi.g, phiI.g, phi0.g, Sr.g, gr, phiI.pinf);
  NumCalc::additivePIntegral<rti>(phi.g, phi0.g, phiI.f, Sr.f, gr, phiI.pinf);
  NumCalc::additivePIntegral<rti>(phi.g, phi0.g, phiI.g, Sr.g, gr, phiI.pinf);
  phi *= (alpha / w2);
  // std::cout << phi * phi << "\n";
}

} // namespace Adams
} // namespace DiracODE
