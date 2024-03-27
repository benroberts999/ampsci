#include "InhomogenousGreens.hpp"
#include "BoundState.hpp"
#include "DiracODE.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <vector>
/*

Solve inhomogenous Dirac equation:
(H_0 + v - e)Fa = S

S (source) is spinor

*/

namespace DiracODE {

DiracSpinor
solve_inhomog(const int kappa, const double en, const std::vector<double> &v,
              const std::vector<double> &H_mag, const double alpha,
              const DiracSpinor &source, const DiracSpinor *const VxFa,
              const DiracSpinor *const Fa0, double zion, double mass) {
  auto Fa = DiracSpinor(0, kappa, source.grid_sptr());
  solve_inhomog(Fa, en, v, H_mag, alpha, source, VxFa, Fa0, zion, mass);
  return Fa;
}

//==============================================================================
void solve_inhomog(DiracSpinor &Fa, const double en,
                   const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source, const DiracSpinor *const VxFa,
                   const DiracSpinor *const Fa0, double zion, double mass)
// NOTE: returns NON-normalised function!
{
  auto Fzero = DiracSpinor(Fa.n(), Fa.kappa(), Fa.grid_sptr());
  auto Finf = DiracSpinor(Fa.n(), Fa.kappa(), Fa.grid_sptr());
  solve_inhomog(Fa, Fzero, Finf, en, v, H_mag, alpha, source, VxFa, Fa0, zion,
                mass);
}
//------------------------------------------------------------------------------
void solve_inhomog(DiracSpinor &Fa, DiracSpinor &Fzero, DiracSpinor &Finf,
                   const double en, const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source, const DiracSpinor *const VxFa,
                   const DiracSpinor *const Fa0, double zion, double mass)
// Overload of the above. Faster, since doesn't need to allocate for Fzero and
// Finf
{
  regularAtOrigin(Fzero, en, v, H_mag, alpha, VxFa, Fa0, zion, mass);
  regularAtInfinity(Finf, en, v, H_mag, alpha, VxFa, Fa0, zion, mass);
  Fa.en() = en;
  Internal::GreenSolution(Fa, Finf, Fzero, alpha, source);
}

namespace Internal {
//==============================================================================
void GreenSolution(DiracSpinor &Fa, const DiracSpinor &Finf,
                   const DiracSpinor &Fzero, const double alpha,
                   const DiracSpinor &Sr) {

  // Wronskian: Should be independent of r
  const auto pp = std::size_t(0.65 * double(Finf.max_pt()));
  auto w2 = (Finf.f(pp) * Fzero.g(pp) - Fzero.f(pp) * Finf.g(pp));
  int f = 1;
  for (auto pt = pp - 20; pt <= pp + 20; ++pt) {
    ++f;
    w2 += (Finf.f(pt) * Fzero.g(pt) - Fzero.f(pt) * Finf.g(pt));
  }
  w2 /= f;

  // std::vector<double> invW2(Finf.f().size());
  // for (std::size_t i = 0; i < Finf.max_pt(); ++i) {
  //   auto w2_i = alpha / (Finf.f(i) * Fzero.g(i) - Fzero.f(i) * Finf.g(i));
  //   auto w2_p = alpha / w2;
  //   invW2[i] = (1.0 * w2_i + 5.0 * w2_p) / 6.0;
  // }

  // save typing:
  const auto &gr = Fa.grid();
  constexpr auto ztr = NumCalc::zero_to_r;
  constexpr auto rti = NumCalc::r_to_inf;

  Fa.max_pt() = gr.num_points();
  Fa *= 0.0;
  Fa.max_pt() = Finf.max_pt();

  NumCalc::additivePIntegral<ztr>(Fa.f(), Finf.f(), Fzero.f(), Sr.f(), gr,
                                  Finf.max_pt());
  NumCalc::additivePIntegral<ztr>(Fa.f(), Finf.f(), Fzero.g(), Sr.g(), gr,
                                  Finf.max_pt());
  NumCalc::additivePIntegral<rti>(Fa.f(), Fzero.f(), Finf.f(), Sr.f(), gr,
                                  Finf.max_pt());
  NumCalc::additivePIntegral<rti>(Fa.f(), Fzero.f(), Finf.g(), Sr.g(), gr,
                                  Finf.max_pt());
  NumCalc::additivePIntegral<ztr>(Fa.g(), Finf.g(), Fzero.f(), Sr.f(), gr,
                                  Finf.max_pt());
  NumCalc::additivePIntegral<ztr>(Fa.g(), Finf.g(), Fzero.g(), Sr.g(), gr,
                                  Finf.max_pt());
  NumCalc::additivePIntegral<rti>(Fa.g(), Fzero.g(), Finf.f(), Sr.f(), gr,
                                  Finf.max_pt());
  NumCalc::additivePIntegral<rti>(Fa.g(), Fzero.g(), Finf.g(), Sr.g(), gr,
                                  Finf.max_pt());
  Fa *= (alpha / w2);

  // const auto R0 =
  //     qip::add(qip::multiply(Fzero.f(), Sr.f), qip::multiply(Fzero.g(),
  //     Sr.g));
  // const auto Ri =
  //     qip::add(qip::multiply(Finf.f(), Sr.f), qip::multiply(Finf.g(), Sr.g));
  // NumCalc::additivePIntegral<ztr>(Fa.f(), Finf.f(), R0, gr, Finf.max_pt());
  // NumCalc::additivePIntegral<rti>(Fa.f(), Fzero.f(), Ri, gr, Finf.max_pt());
  // NumCalc::additivePIntegral<ztr>(Fa.g(), Finf.g(), R0, gr, Finf.max_pt());
  // NumCalc::additivePIntegral<rti>(Fa.g(), Fzero.g(), Ri, gr, Finf.max_pt());
  // Fa *= (alpha / w2);
}

} // namespace Internal
} // namespace DiracODE
