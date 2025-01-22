#include "InhomogenousGreens.hpp"
#include "BoundState.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "include.hpp"
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
void GreenSolution(DiracSpinor &Fa, const DiracSpinor &wi,
                   const DiracSpinor &w0, const double alpha,
                   const DiracSpinor &s) {

  /*
      Calculates:
      Fa(r) = (1/w) * ( ∫_0^r w0(x)s(x)dx + ∫_r^∞ wi(x)s(x)dx )
  */

  // Wronskian: Should be independent of r
  const auto pp = std::size_t(0.65 * double(wi.max_pt()));
  auto w2 = (wi.f(pp) * w0.g(pp) - w0.f(pp) * wi.g(pp));
  int count = 1;
  for (auto pt = pp - 20; pt <= pp + 20; ++pt) {
    ++count;
    w2 += (wi.f(pt) * w0.g(pt) - w0.f(pt) * wi.g(pt));
  }
  w2 /= count;

  // R-dependent wronskian
  const auto wr = [&](std::size_t i) {
    const auto tmp = 0.5 * (w2 + wi.f(i) * w0.g(i) - w0.f(i) * wi.g(i));
    return tmp == 0.0 ? w2 : tmp;
  };

  // save typing:
  const auto &gr = Fa.grid();
  const auto irmax = std::max(wi.max_pt(), Fa.max_pt());

  // clear existing solution
  // (don't need to clear between [0,rmax], it is overwritten below)
  Fa.min_pt() = 0;
  Fa.max_pt() = irmax;
  Fa.zero_boundaries();

  const auto du = gr.du();
  const auto num_points = gr.num_points();

  // Quadrature integration weights:
  const auto dr = [&](std::size_t i) {
    if (i < NumCalc::Nquad) {
      return NumCalc::dq_inv * NumCalc::cq[i] * gr.drdu(i);
    }
    if (i < num_points - NumCalc::Nquad)
      return gr.drdu(i);
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1] * gr.drdu(i);
  };

  double A = 0.5 * (w0.f(0) * s.f(0) + w0.g(0) * s.g(0)) * dr(0);
  double B = 0.0;

  Fa.f(0) = (wi.f(0) * A);
  Fa.g(0) = (wi.g(0) * A);

  for (std::size_t i = 1; i < irmax; ++i) {
    A += 0.5 *
         ((w0.f(i - 1) * s.f(i - 1) + w0.g(i - 1) * s.g(i - 1)) * dr(i - 1) +
          (w0.f(i) * s.f(i) + w0.g(i) * s.g(i)) * dr(i));

    Fa.f(i) = wi.f(i) * A;
    Fa.g(i) = wi.g(i) * A;
  }

  for (std::size_t ii = irmax - 1; ii >= 1; --ii) {
    const auto i = ii - 1;

    B += 0.5 *
         ((wi.f(i + 1) * s.f(i + 1) + wi.g(i + 1) * s.g(i + 1)) * dr(i + 1) +
          (wi.f(i) * s.f(i) + wi.g(i) * s.g(i)) * dr(i));

    Fa.f(i) += w0.f(i) * B;
    Fa.g(i) += w0.g(i) * B;
  }

  // Wronskian - Should be r=independant, but allow to be more general
  for (std::size_t i = 0; i < irmax; ++i) {
    Fa.f(i) *= alpha * du / wr(i);
    Fa.g(i) *= alpha * du / wr(i);
  }
}

} // namespace Internal
} // namespace DiracODE
