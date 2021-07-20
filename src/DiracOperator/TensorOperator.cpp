#include "DiracOperator/TensorOperator.hpp"
#include "Angular/Angular_369j.hpp"
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

  const auto &df = (diff_order == 0) ? Fb.f()
                                     : NumCalc::derivative(Fb.f(), gr.drdu(),
                                                           gr.du(), diff_order);
  const auto &dg = (diff_order == 0) ? Fb.g()
                                     : NumCalc::derivative(Fb.g(), gr.drdu(),
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

  return Fa * radial_rhs(Fa.k, Fb);

  // const auto &gr = *(Fb.rgrid);
  // if (isZero(Fa.k, Fb.k))
  //   return 0.0;
  //
  // const auto cff = angularCff(Fa.k, Fb.k);
  // const auto cgg = angularCgg(Fa.k, Fb.k);
  // const auto cfg = angularCfg(Fa.k, Fb.k);
  // const auto cgf = angularCgf(Fa.k, Fb.k);
  //
  // auto p0 = std::max(Fa.min_pt(), Fb.min_pt());
  // auto pinf = std::min(Fa.max_pt(), Fb.max_pt());
  // double radint = 0.0;
  //
  // // df is either just fb, or dFb/dr
  // const auto &df = (diff_order == 0)
  //                      ? Fb.f
  //                      : NumCalc::derivative(Fb.f(), gr.drdu(), gr.du(),
  //                      diff_order);
  // const auto &dg = (diff_order == 0)
  //                      ? Fb.g
  //                      : NumCalc::derivative(Fb.g(), gr.drdu(), gr.du(),
  //                      diff_order);
  //
  // if (m_vec.empty()) {
  //   if (cff != 0.0)
  //     radint += NumCalc::integrate(cff, p0, pinf, Fa.f(), df, gr.drdu());
  //   if (cfg != 0.0)
  //     radint += NumCalc::integrate(cfg, p0, pinf, Fa.f(), dg, gr.drdu());
  //   if (cgf != 0.0)
  //     radint += NumCalc::integrate(cgf, p0, pinf, Fa.g(), df, gr.drdu());
  //   if (cgg != 0.0)
  //     radint += NumCalc::integrate(cgg, p0, pinf, Fa.g(), dg, gr.drdu());
  // } else {
  //   if (cff != 0.0)
  //     radint += NumCalc::integrate(cff, p0, pinf, Fa.f(), df, gr.drdu(),
  //     m_vec);
  //   if (cfg != 0.0)
  //     radint += NumCalc::integrate(cfg, p0, pinf, Fa.f(), dg, gr.drdu(),
  //     m_vec);
  //   if (cgf != 0.0)
  //     radint += NumCalc::integrate(cgf, p0, pinf, Fa.g(), df, gr.drdu(),
  //     m_vec);
  //   if (cgg != 0.0)
  //     radint += NumCalc::integrate(cgg, p0, pinf, Fa.g(), dg, gr.drdu(),
  //     m_vec);
  // }
  //
  // return m_constant * radint * gr.du();
}

//******************************************************************************
//******************************************************************************
//******************************************************************************
IntM4x4 IntM4x4::operator*(const IntM4x4 &other) const {
  bool imag = false;
  int sign = 1;
  if (this->imaginary && other.imaginary) {
    sign = -1;
  } else if (this->imaginary || other.imaginary) {
    imag = true;
  }
  int a = sign * (this->e00 * other.e00 + this->e01 * other.e10);
  int b = sign * (this->e00 * other.e01 + this->e01 * other.e11);
  int c = sign * (this->e10 * other.e00 + this->e11 * other.e10);
  int d = sign * (this->e10 * other.e01 + this->e11 * other.e11);
  return IntM4x4(a, b, c, d, imag);
}
IntM4x4 IntM4x4::operator+(const IntM4x4 &other) const {
  bool imag = false;
  if (this->imaginary && other.imaginary) {
    imag = true;
  } else if (this->imaginary || other.imaginary) {
    std::cerr << "FAIL 50 in ScalarOperator_old. Cannot mix real and imaginary "
                 "Dirac Matrices! \n";
    std::abort();
  }
  int a = this->e00 + other.e00;
  int b = this->e01 + other.e01;
  int c = this->e10 + other.e10;
  int d = this->e11 + other.e11;
  return IntM4x4(a, b, c, d, imag);
}
IntM4x4 IntM4x4::operator-(const IntM4x4 &other) const {
  bool imag = false;
  if (this->imaginary && other.imaginary) {
    imag = true;
  } else if (this->imaginary || other.imaginary) {
    std::cerr << "FAIL 50 in ScalarOperator_old. Cannot mix real and imaginary "
                 "Dirac Matrices! \n";
    std::abort();
  }
  int a = this->e00 - other.e00;
  int b = this->e01 - other.e01;
  int c = this->e10 - other.e10;
  int d = this->e11 - other.e11;
  return IntM4x4(a, b, c, d, imag);
}

} // namespace DiracOperator
