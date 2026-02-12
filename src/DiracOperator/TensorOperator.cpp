#include "DiracOperator/TensorOperator.hpp"
#include "Angular/Wigner369j.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace DiracOperator {

//==============================================================================
bool TensorOperator::isZero(const int ka, int kb) const {
  // checks m_rank and m_parity

  // if (m_rank < std::abs(Angular::twoj_k(ka) - Angular::twoj_k(kb)) / 2)
  //   return true;

  if (Angular::triangle(Angular::twoj_k(ka), Angular::twoj_k(kb), 2 * m_rank) ==
      0)
    return true;

  if ((m_parity == Parity::even) !=
      (Angular::parity_k(ka) == Angular::parity_k(kb)))
    return true;

  return false;
  // can remove this??
  // return (angularF(ka, kb) == 0);
}

bool TensorOperator::isZero(const DiracSpinor &Fa,
                            const DiracSpinor &Fb) const {
  return isZero(Fa.kappa(), Fb.kappa());
}

//==============================================================================
double TensorOperator::rme3js(const int twoja, const int twojb, int two_mb,
                              int two_q) const {
  // rme3js = (-1)^{ja-ma} (ja, k, jb,\ -ma, q, mb)
  const auto two_ma = two_mb + two_q; // -ma + mb + q = 0;
  // sign = (-1)^(ja - ma)
  const auto sign = ((twoja - two_ma) / 2) % 2 == 0 ? 1 : -1;
  return sign *
         Angular::threej_2(twoja, 2 * m_rank, twojb, -two_ma, two_q, two_mb);
}

DiracSpinor TensorOperator::reduced_rhs(const int ka,
                                        const DiracSpinor &Fb) const {
  return angularF(ka, Fb.kappa()) * radial_rhs(ka, Fb);
}

DiracSpinor TensorOperator::reduced_lhs(const int ka,
                                        const DiracSpinor &Fb) const {
  const int s = imaginaryQ() ? -1 : 1;
  const auto x = Angular::evenQ_2(Angular::twoj_k(ka) - Fb.twoj()) ? s : -s;
  return (x * angularF(ka, Fb.kappa())) * radial_rhs(ka, Fb);
}

double TensorOperator::reducedME(const DiracSpinor &Fa,
                                 const DiracSpinor &Fb) const {
  return angularF(Fa.kappa(), Fb.kappa()) * radialIntegral(Fa, Fb);
}

double TensorOperator::fullME(const DiracSpinor &Fa, const DiracSpinor &Fb,
                              std::optional<int> two_ma,
                              std::optional<int> two_mb,
                              std::optional<int> two_q) const {

  const auto tma = two_ma ? *two_ma : std::min(Fa.twoj(), Fb.twoj());
  const auto tmb = two_mb ? *two_mb : tma;
  const auto tqq = two_q ? *two_q : 0;

  const auto sign = Angular::neg1pow_2(Fa.twoj() - tma);
  const auto factor = sign * Angular::threej_2(Fa.twoj(), 2 * m_rank, Fb.twoj(),
                                               -tma, tqq, tmb);

  return factor * reducedME(Fa, Fb);
}

//==============================================================================
void TensorOperator::scale(double lambda) { m_constant *= lambda; }

//==============================================================================

DiracSpinor TensorOperator::radial_rhs(const int kappa_a,
                                       const DiracSpinor &Fb) const {
  // Fa * radial_rhs(Fa.kappa(),Fb) = h.radialIntegral(Fa, Fb)

  const auto &gr = Fb.grid();
  DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
  dF.min_pt() = Fb.min_pt();
  dF.max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.kappa())) {
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.min_pt();
    return dF;
  }

  const auto &df =
      (m_diff_order == 0) ?
          Fb.f() :
          NumCalc::derivative(Fb.f(), gr.drdu(), gr.du(), m_diff_order);
  const auto &dg =
      (m_diff_order == 0) ?
          Fb.g() :
          NumCalc::derivative(Fb.g(), gr.drdu(), gr.du(), m_diff_order);

  const auto cff = angularCff(kappa_a, Fb.kappa());
  const auto cgg = angularCgg(kappa_a, Fb.kappa());
  const auto cfg = angularCfg(kappa_a, Fb.kappa());
  const auto cgf = angularCgf(kappa_a, Fb.kappa());

  if (m_vec.empty()) {
    for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
      dF.f(i) = m_constant * (cff * df[i] + cfg * dg[i]);
      dF.g(i) = m_constant * (cgf * df[i] + cgg * dg[i]);
    }
  } else {
    for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
      dF.f(i) = m_constant * m_vec[i] * (cff * df[i] + cfg * dg[i]);
      dF.g(i) = m_constant * m_vec[i] * (cgf * df[i] + cgg * dg[i]);
    }
  }

  return dF;
}

//==============================================================================
double TensorOperator::radialIntegral(const DiracSpinor &Fa,
                                      const DiracSpinor &Fb) const {

  const int kappa_a = Fa.kappa();
  if (isZero(kappa_a, Fb.kappa())) {
    return 0.0;
  }

  const auto &gr = Fb.grid();

  const auto p0 = std::max(Fa.min_pt(), Fb.min_pt());
  const auto pf = std::min(Fa.max_pt(), Fb.max_pt());

  const auto &df =
      (m_diff_order == 0) ?
          Fb.f() :
          NumCalc::derivative(Fb.f(), gr.drdu(), gr.du(), m_diff_order);
  const auto &dg =
      (m_diff_order == 0) ?
          Fb.g() :
          NumCalc::derivative(Fb.g(), gr.drdu(), gr.du(), m_diff_order);

  const auto cff = angularCff(kappa_a, Fb.kappa());
  const auto cgg = angularCgg(kappa_a, Fb.kappa());
  const auto cfg = angularCfg(kappa_a, Fb.kappa());
  const auto cgf = angularCgf(kappa_a, Fb.kappa());

  const auto Rff =
      cff == 0.0 ?
          0.0 :
      m_vec.empty() ?
          cff * NumCalc::integrate(1.0, p0, pf, Fa.f(), df, gr.drdu()) :
          cff * NumCalc::integrate(1.0, p0, pf, m_vec, Fa.f(), df, gr.drdu());
  const auto Rfg =
      cfg == 0.0 ?
          0.0 :
      m_vec.empty() ?
          cfg * NumCalc::integrate(1.0, p0, pf, Fa.f(), dg, gr.drdu()) :
          cfg * NumCalc::integrate(1.0, p0, pf, m_vec, Fa.f(), dg, gr.drdu());
  const auto Rgf =
      cgf == 0.0 ?
          0.0 :
      m_vec.empty() ?
          cgf * NumCalc::integrate(1.0, p0, pf, Fa.g(), df, gr.drdu()) :
          cgf * NumCalc::integrate(1.0, p0, pf, m_vec, Fa.g(), df, gr.drdu());
  const auto Rgg =
      cgg == 0.0 ?
          0.0 :
      m_vec.empty() ?
          cgg * NumCalc::integrate(1.0, p0, pf, Fa.g(), dg, gr.drdu()) :
          cgg * NumCalc::integrate(1.0, p0, pf, m_vec, Fa.g(), dg, gr.drdu());

  return m_constant * gr.du() * (Rff + Rfg + Rgf + Rgg);
}

//******************************************************************************
// Helper functions: Useful for several operators
//******************************************************************************

// Pab function: Int[ (fa*gb + pm*ga*fb) * t(r) , dr]. pm = +/-1 (usually)
double Pab(double pm, const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb) {
  if (pm == -1 && &Fa == &Fb)
    return 0.0;
  const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
  const auto pf = std::min(Fa.max_pt(), Fb.max_pt());
  const auto &drdu = Fb.grid().drdu();
  const auto fg = NumCalc::integrate(1.0, pi, pf, t, Fa.f(), Fb.g(), drdu);
  const auto gf = &Fa == &Fb ?
                      fg :
                      NumCalc::integrate(1.0, pi, pf, t, Fa.g(), Fb.f(), drdu);
  return (fg + pm * gf) * Fb.grid().du();
}

// Rab function: Int[ (fa*fb + pm*ga*gb) * t(r) , dr]. pm = +/-1 (usually)
double Rab(double pm, const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb) {
  const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
  const auto pf = std::min(Fa.max_pt(), Fb.max_pt());
  const auto &drdu = Fb.grid().drdu();
  const auto ff = NumCalc::integrate(1.0, pi, pf, t, Fa.f(), Fb.f(), drdu);
  const auto gg = NumCalc::integrate(1.0, pi, pf, t, Fa.g(), Fb.g(), drdu);
  return (ff + pm * gg) * Fb.grid().du();
}

// Pab_rhs function: dF_ab += t(r) * (g, pm*f) - note, uses +=, so can combine.
// Ensure empty to begin.
void Pab_rhs(double pm, const std::vector<double> &t, DiracSpinor *dF,
             const DiracSpinor &Fb, double a) {
  for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
    dF->f(i) += a * t[i] * Fb.g(i);
    dF->g(i) += a * pm * t[i] * Fb.f(i);
  }
}

// Rab_rhs function: dF_ab += t(r) * (f, pm*g) - note, uses +=, so can combine.
// Ensure empty to begin.
void Rab_rhs(double pm, const std::vector<double> &t, DiracSpinor *dF,
             const DiracSpinor &Fb, double a) {
  for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
    dF->f(i) += a * t[i] * Fb.f(i);
    dF->g(i) += a * pm * t[i] * Fb.g(i);
  }
}

// Vab function: Int[ (fa*gb ) * t(r) , dr].
double Vab(const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb) {
  const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
  const auto pf = std::min(Fa.max_pt(), Fb.max_pt());
  const auto &drdu = Fb.grid().drdu();
  const auto fg = NumCalc::integrate(1.0, pi, pf, t, Fa.f(), Fb.g(), drdu);
  return fg * Fb.grid().du();
}

// Wab function: Int[ (ga*fb ) * t(r) , dr].
double Wab(const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb) {
  const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
  const auto pf = std::min(Fa.max_pt(), Fb.max_pt());
  const auto &drdu = Fb.grid().drdu();
  const auto fg = NumCalc::integrate(1.0, pi, pf, t, Fa.g(), Fb.f(), drdu);
  return fg * Fb.grid().du();
}

double Gab(const std::vector<double> &t, const DiracSpinor &Fa,
           const DiracSpinor &Fb) {
  const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
  const auto pf = std::min(Fa.max_pt(), Fb.max_pt());
  const auto &drdu = Fb.grid().drdu();
  const auto gg = NumCalc::integrate(1.0, pi, pf, t, Fa.g(), Fb.g(), drdu);
  return gg * Fb.grid().du();
}

void Gab_rhs(const std::vector<double> &t, DiracSpinor *dF,
             const DiracSpinor &Fb, double a) {
  for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
    // dF->f(i) += a * t[i] * Fb.f(i);
    dF->g(i) += a * t[i] * Fb.g(i);
  }
}

//******************************************************************************
double Gab(const DiracSpinor &Fa, const DiracSpinor &Fb) {
  const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
  const auto pf = std::min(Fa.max_pt(), Fb.max_pt());
  const auto &drdu = Fb.grid().drdu();
  const auto gg = NumCalc::integrate(1.0, pi, pf, Fa.g(), Fb.g(), drdu);
  return gg * Fb.grid().du();
}

// Gab_rhs function (constant version): dF_ab += g - note, uses +=, so can combine.
// Ensure empty to begin.
void Gab_rhs(DiracSpinor *dF, const DiracSpinor &Fb, double a) {
  for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
    dF->g(i) += a * Fb.g(i);
  }
}

//******************************************************************************
// Versions for constant t(r) = c

// Pab function: Int[ (fa*gb + pm*ga*fb) , dr]. pm = +/-1 (usually)
double Pab(double pm, const DiracSpinor &Fa, const DiracSpinor &Fb) {
  if (pm == -1 && &Fa == &Fb)
    return 0.0;
  const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
  const auto pf = std::min(Fa.max_pt(), Fb.max_pt());
  const auto &drdu = Fb.grid().drdu();
  const auto fg = NumCalc::integrate(1.0, pi, pf, Fa.f(), Fb.g(), drdu);
  const auto gf =
      &Fa == &Fb ? fg : NumCalc::integrate(1.0, pi, pf, Fa.g(), Fb.f(), drdu);
  return (fg + pm * gf) * Fb.grid().du();
}

double Rab(double pm, const DiracSpinor &Fa, const DiracSpinor &Fb) {

  const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
  const auto pf = std::min(Fa.max_pt(), Fb.max_pt());
  const auto &drdu = Fb.grid().drdu();

  // NO! That;s what G is for!
  // ONLY if kappa_a = kappa_b!!!
  // // Use orthogonality of Fa and Fb:
  // // (ff + a*gg) = (ff + gg + [a-1]*gg) = [a-1]*gg
  // if (pm == 1.0)
  // return 0.0;
  // const auto gg = NumCalc::integrate(1.0, pi, pf, Fa.g(), Fb.g(), drdu);
  // return (pm - 1.0) * gg * Fb.grid().du();

  const auto ff = NumCalc::integrate(1.0, pi, pf, Fa.f(), Fb.f(), drdu);
  const auto gg = NumCalc::integrate(1.0, pi, pf, Fa.g(), Fb.g(), drdu);
  return (ff + pm * gg) * Fb.grid().du();
}

// Pab_rhs function: dF_ab += t(r) * (g, pm*f) - note, uses +=, so can combine.
// Ensure empty to begin.
void Pab_rhs(double pm, DiracSpinor *dF, const DiracSpinor &Fb, double a) {
  for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
    dF->f(i) += a * Fb.g(i);
    dF->g(i) += a * pm * Fb.f(i);
  }
}

// Rab_rhs function: dF_ab += t(r) * (f, pm*g) - note, uses +=, so can combine.
// Ensure empty to begin.
void Rab_rhs(double pm, DiracSpinor *dF, const DiracSpinor &Fb, double a) {

  // NO! That;s what G is for!
  // ONLY if kappa_a = kappa_b!!!
  // Use orthogonality of Fa and Fb:
  // (ff + a*gg) = (ff + gg + [a-1]*gg) = [a-1]*gg
  // for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
  //   dF->f(i) += 0.0;
  //   dF->g(i) += a * (pm - 1.0) * Fb.g(i);
  // }

  for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
    dF->f(i) += a * Fb.f(i);
    dF->g(i) += a * pm * Fb.g(i);
  }
}

} // namespace DiracOperator
