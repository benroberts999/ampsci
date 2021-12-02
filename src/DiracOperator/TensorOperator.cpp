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

//******************************************************************************
bool TensorOperator::isZero(const int ka, int kb) const {
  // checks m_rank and m_parity
  if (m_rank < std::abs(Angular::twoj_k(ka) - Angular::twoj_k(kb)) / 2)
    return true;
  if ((m_parity == Parity::even) !=
      (Angular::parity_k(ka) == Angular::parity_k(kb)))
    return true;
  if (angularF(ka, kb) == 0)
    return true;
  return false; /*may still be zero*/
}
bool TensorOperator::isZero(const DiracSpinor &Fa,
                            const DiracSpinor &Fb) const {
  return isZero(Fa.k, Fb.k);
}

std::string TensorOperator::rme_symbol(const DiracSpinor &Fa,
                                       const DiracSpinor &Fb) const {
  return std::string("<") + Fa.shortSymbol() + "||h||" + Fb.shortSymbol() + ">";
}
std::string TensorOperator::R_symbol(const DiracSpinor &Fa,
                                     const DiracSpinor &Fb) const {
  return std::string("R(") + Fa.shortSymbol() + "," + Fb.shortSymbol() + ")";
}

//******************************************************************************
double TensorOperator::rme3js(const int twoja, const int twojb, int two_mb,
                              int two_q) const {
  // rme3js = (-1)^{ja-ma} (ja, k, jb,\ -ma, q, mb)
  const auto two_ma = two_mb - two_q; // -ma + mb + q = 0;
  // sig = (-1)^(ja - ma)
  const auto sig = ((twoja - two_ma) / 2) % 2 == 0 ? 1 : -1;
  return sig *
         Angular::threej_2(twoja, 2 * m_rank, twojb, -two_ma, two_q, two_mb);
}

DiracSpinor TensorOperator::reduced_rhs(const int ka,
                                        const DiracSpinor &Fb) const {
  return angularF(ka, Fb.k) * radial_rhs(ka, Fb);
}

DiracSpinor TensorOperator::reduced_lhs(const int ka,
                                        const DiracSpinor &Fb) const {
  const int s = imaginaryQ() ? -1 : 1;
  const auto x = Angular::evenQ_2(Angular::twoj_k(ka) - Fb.twoj()) ? s : -s;
  return (x * angularF(ka, Fb.k)) * radial_rhs(ka, Fb);
}

double TensorOperator::reducedME(const DiracSpinor &Fa,
                                 const DiracSpinor &Fb) const {
  return angularF(Fa.k, Fb.k) * radialIntegral(Fa, Fb);
}

//******************************************************************************

DiracSpinor TensorOperator::radial_rhs(const int kappa_a,
                                       const DiracSpinor &Fb) const {
  // Fa * radial_rhs(Fa.k,Fb) = h.radialIntegral(Fa, Fb)

  const auto &gr = *(Fb.rgrid);
  DiracSpinor dF(0, kappa_a, Fb.rgrid);
  dF.set_min_pt() = Fb.min_pt();
  dF.set_max_pt() = Fb.max_pt();

  if (isZero(kappa_a, Fb.k)) {
    dF.set_min_pt() = Fb.min_pt();
    dF.set_max_pt() = Fb.min_pt();
    return dF;
  }

  const auto &df = (diff_order == 0) ? Fb.f() :
                                       NumCalc::derivative(Fb.f(), gr.drdu(),
                                                           gr.du(), diff_order);
  const auto &dg = (diff_order == 0) ? Fb.g() :
                                       NumCalc::derivative(Fb.g(), gr.drdu(),
                                                           gr.du(), diff_order);

  const auto cff = angularCff(kappa_a, Fb.k);
  const auto cgg = angularCgg(kappa_a, Fb.k);
  const auto cfg = angularCfg(kappa_a, Fb.k);
  const auto cgf = angularCgf(kappa_a, Fb.k);

  if (m_vec.empty()) {
    for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
      dF.set_f(i) = m_constant * (cff * df[i] + cfg * dg[i]);
      dF.set_g(i) = m_constant * (cgf * df[i] + cgg * dg[i]);
    }
  } else {
    for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
      dF.set_f(i) = m_constant * m_vec[i] * (cff * df[i] + cfg * dg[i]);
      dF.set_g(i) = m_constant * m_vec[i] * (cgf * df[i] + cgg * dg[i]);
    }
  }

  return dF;
}

//******************************************************************************
double TensorOperator::radialIntegral(const DiracSpinor &Fa,
                                      const DiracSpinor &Fb) const {

  // nb: faster not to do this, but nicer this way
  return Fa * radial_rhs(Fa.k, Fb);
}

} // namespace DiracOperator
