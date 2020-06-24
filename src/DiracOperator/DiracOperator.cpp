#include "DiracOperator/DiracOperator.hpp"
#include "Angular/Angular_369j.hpp"
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
                              int two_q) const
// rme3js = (-1)^{ja-ma} (ja, k, jb,\ -ma, q, mb)
{
  auto two_ma = two_mb - two_q; // -ma + mb + q = 0;
  // sig = (-1)^(ja - ma)
  auto sig = ((twoja - two_ma) / 2) % 2 == 0 ? 1 : -1;
  return sig *
         Angular::threej_2(twoja, 2 * m_rank, twojb, -two_ma, two_q, two_mb);
}

DiracSpinor TensorOperator::reduced_rhs(const int ka,
                                        const DiracSpinor &Fb) const {
  return angularF(ka, Fb.k) * radial_rhs(ka, Fb);
}

DiracSpinor TensorOperator::reduced_lhs(const int ka,
                                        const DiracSpinor &Fb) const {
  int s = imaginaryQ() ? -1 : 1;
  auto x = Angular::evenQ_2(Angular::twoj_k(ka) - Fb.twoj()) ? s : -s;
  return (x * angularF(ka, Fb.k)) * radial_rhs(ka, Fb);
}

double TensorOperator::radialIntegral(const DiracSpinor &Fa,
                                      const DiracSpinor &Fb) const {
  if (isZero(Fa.k, Fb.k))
    return 0.0;
  return (Fa * radial_rhs(Fa.k, Fb));
}
double TensorOperator::reducedME(const DiracSpinor &Fa,
                                 const DiracSpinor &Fb) const {
  return angularF(Fa.k, Fb.k) * radialIntegral(Fa, Fb);
}

//******************************************************************************
DiracSpinor TensorOperator::radial_rhs(const int kappa_a,
                                       const DiracSpinor &Fb) const
// psi1 * dPsi2 = h.radialIntegral(psi1,psi2)
// Because of angular factor, _may_ depend on kappa of 'lhs'
// Usually, will not over-write this. But in some cases, might be better
// (eg, M1)
{

  const auto &gr = *(Fb.rgrid);
  DiracSpinor dPsi(0, kappa_a, Fb.rgrid);
  if (isZero(kappa_a, Fb.k))
    return dPsi;

  // Strangeness here to account for possible derivatives
  // (I am trying to avoid doing a copy when no derivative)
  // Copy is unavoidable when calcing derivative, but that's rare!
  const std::vector<double> *rhs_f = &(Fb.f);
  const std::vector<double> *rhs_g = &(Fb.g);
  std::unique_ptr<const std::vector<double>> dummy_df = nullptr;
  std::unique_ptr<const std::vector<double>> dummy_dg = nullptr;
  if (diff_order > 0) {
    // rhs_f is either a pointer to F_input, OR (in case of deriv, F')
    dummy_df = std::make_unique<std::vector<double>>(
        NumCalc::derivative(Fb.f, gr.drdu, gr.du, diff_order));
    dummy_dg = std::make_unique<std::vector<double>>(
        NumCalc::derivative(Fb.g, gr.drdu, gr.du, diff_order));
    rhs_f = dummy_df.get();
    rhs_g = dummy_dg.get();
  }

  const auto cff = angularCff(kappa_a, Fb.k);
  const auto cgg = angularCgg(kappa_a, Fb.k);
  const auto cfg = angularCfg(kappa_a, Fb.k);
  const auto cgf = angularCgf(kappa_a, Fb.k);
  for (unsigned i = 0; i < Fb.pinf; i++) {
    dPsi.f[i] = m_constant * (cff * (*rhs_f)[i] + cfg * (*rhs_g)[i]);
    dPsi.g[i] = m_constant * (cgf * (*rhs_f)[i] + cgg * (*rhs_g)[i]);
  }
  if (!m_vec.empty()) {
    dPsi *= m_vec;
  }

  return dPsi;
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
