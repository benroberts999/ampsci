#include "Adams_Greens.hpp"
#include "Adams_bound.hpp"
#include "DiracODE.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>
/*

Solve inhomogenous Dirac equation:
(H_0 + v - e)phi = S

S (source) is spinor

*/

namespace DiracODE {

DiracSpinor solve_inhomog(const int kappa, const double en,
                          const std::vector<double> &v,
                          const std::vector<double> &H_mag, const double alpha,
                          const DiracSpinor &source) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__, "0");
  auto phi = DiracSpinor(0, kappa, source.rgrid);
  solve_inhomog(phi, en, v, H_mag, alpha, source);
  return phi;
}

//******************************************************************************
void solve_inhomog(DiracSpinor &phi, const double en,
                   const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source)
// NOTE: returns NON-normalised function!
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__, "a");
  auto phi0 = DiracSpinor(phi.n, phi.k, phi.rgrid);
  auto phiI = DiracSpinor(phi.n, phi.k, phi.rgrid);
  solve_inhomog(phi, phi0, phiI, en, v, H_mag, alpha, source);
}
//------------------------------------------------------------------------------
void solve_inhomog(DiracSpinor &phi, DiracSpinor &phi0, DiracSpinor &phiI,
                   const double en, const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source)
// Overload of the above. Faster, since doesn't need to allocate for phi0 and
// phiI
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__, "b");
  regularAtOrigin(phi0, en, v, H_mag, alpha);
  regularAtInfinity(phiI, en, v, H_mag, alpha);
  Adams::GreenSolution(phi, phiI, phi0, alpha, source);
}

namespace Adams {
//******************************************************************************
void GreenSolution(DiracSpinor &phi, const DiracSpinor &phiI,
                   const DiracSpinor &phi0, const double alpha,
                   const DiracSpinor &Sr) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // Wronskian:
  auto pp = std::size_t(0.65 * double(phiI.pinf));
  auto w2 = (phiI.f[pp] * phi0.g[pp] - phi0.f[pp] * phiI.g[pp]);

  // std::vector<double> invW2(phiI.f.size());
  // for (std::size_t i = 0; i < phiI.pinf; ++i) {
  //   auto w2_i = alpha / (phiI.f[i] * phi0.g[i] - phi0.f[i] * phiI.g[i]);
  //   auto w2_p = alpha / w2;
  //   invW2[i] = (1.0 * w2_i + 5.0 * w2_p) / 6.0;
  // }

  // save typing:
  const auto &gr = *phi.rgrid;
  constexpr auto ztr = NumCalc::zero_to_r;
  constexpr auto rti = NumCalc::r_to_inf;

  phi.pinf = gr.num_points;
  phi *= 0.0;
  phi.pinf = phiI.pinf;
  NumCalc::additivePIntegral<ztr>(phi.f, phiI.f, phi0.f, Sr.f, gr, phiI.pinf);
  NumCalc::additivePIntegral<ztr>(phi.f, phiI.f, phi0.g, Sr.g, gr, phiI.pinf);
  NumCalc::additivePIntegral<rti>(phi.f, phi0.f, phiI.f, Sr.f, gr, phiI.pinf);
  NumCalc::additivePIntegral<rti>(phi.f, phi0.f, phiI.g, Sr.g, gr, phiI.pinf);
  NumCalc::additivePIntegral<ztr>(phi.g, phiI.g, phi0.f, Sr.f, gr, phiI.pinf);
  NumCalc::additivePIntegral<ztr>(phi.g, phiI.g, phi0.g, Sr.g, gr, phiI.pinf);
  NumCalc::additivePIntegral<rti>(phi.g, phi0.g, phiI.f, Sr.f, gr, phiI.pinf);
  NumCalc::additivePIntegral<rti>(phi.g, phi0.g, phiI.g, Sr.g, gr, phiI.pinf);
  phi *= (alpha / w2);
  // phi *= invW2;
}

} // namespace Adams
} // namespace DiracODE
