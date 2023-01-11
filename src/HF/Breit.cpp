#include "HF/Breit.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "Maths/Grid.hpp"
#include <iostream>
#include <memory>

namespace HF {

//==============================================================================
DiracSpinor Breit::VbrFa(const DiracSpinor &Fa,
                         const std::vector<DiracSpinor> &core) const {
  DiracSpinor BFa(Fa.n(), Fa.kappa(), Fa.grid_sptr());
  for (const auto &Fb : core) {
    const auto kmin = std::abs(Fb.twoj() - Fa.twoj()) / 2;
    const auto kmax = (Fb.twoj() + Fa.twoj()) / 2;
    for (int k = kmin; k <= kmax; ++k) {
      const auto s = Angular::neg1pow(k);
      if (s == 1)
        BFa += Bkv_bcd(k, Fa.kappa(), Fb, Fb, Fa);
      else
        BFa -= Bkv_bcd(k, Fa.kappa(), Fb, Fb, Fa);
    }
  }
  BFa *= (-1.0 / Fa.twojp1());
  return BFa;
}

//==============================================================================
double Breit::Bk_abcd(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                      const DiracSpinor &Fc, const DiracSpinor &Fd) const {
  return Fa * Bkv_bcd(k, Fa.kappa(), Fb, Fc, Fd);
}
double Breit::BPk_abcd(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                       const DiracSpinor &Fc, const DiracSpinor &Fd) const {
  return Fa * BPkv_bcd(k, Fa.kappa(), Fb, Fc, Fd);
}
double Breit::BWk_abcd(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                       const DiracSpinor &Fc, const DiracSpinor &Fd) const {
  return Fa * BWkv_bcd(k, Fa.kappa(), Fb, Fc, Fd);
}

//==============================================================================
DiracSpinor Breit::Bkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                           const DiracSpinor &Fc, const DiracSpinor &Fd) const {

  DiracSpinor out(0, kappa_v, Fc.grid_sptr());
  out.min_pt() = Fc.min_pt();
  out.max_pt() = Fc.max_pt();

  const auto ka = kappa_v;
  const auto kb = Fb.kappa();
  const auto kc = Fc.kappa();
  const auto kd = Fd.kappa();

  const auto Ckac = Angular::Ck_kk(k, ka, kc);
  const auto Ckbd = Angular::Ck_kk(k, kb, kd);
  const auto Dkac = Angular::Ck_kk(k, -ka, kc);
  const auto Dkbd = Angular::Ck_kk(k, -kb, kd);
  const auto tja = Angular::twoj_k(ka);
  const auto sign = Angular::neg1pow_2(2 * k + tja - Fb.twoj());

  // k=0 in MOP part??
  const auto have_mop = Ckac != 0.0 && Ckbd != 0.0;
  const auto have_n = (Dkac != 0.0) && (Dkbd != 0.0) && (ka + kc != 0) &&
                      (kb + kd != 0) && (k != 0.0);

  // nb: never have both MOP _and_ N ! different parity rule for each k!
  assert(!(have_mop && have_n));

  if (have_mop) {
    // Angular factors
    const auto d_ac = ka - kc;
    const auto d_bd = kb - kd;
    const auto eta_k = k == 0.0 ? 0.0 : d_bd / double(k);
    const auto eta_kp1 = d_bd / double(k + 1);
    const auto eacp1 = d_ac / double(k + 1);
    const auto eac = k == 0.0 ? 0.0 : d_ac / double(k);

    // Calculate the Breit radial screening integrals
    const auto gbk = Breit_gb::single_k_mop(Fb, Fd, k);

    // coeficients (including scaling factors)
    const auto factor = sign * m_scale * Ckac * Ckbd;
    const auto c_m1 = m_M * (k + 1) / double(2 * k + 3);
    const auto c_m2 = m_M * k / double(2 * k - 1);
    const auto c_o1 =
        -m_O * (k + 1) * (k + 1) / double((2 * k + 1) * (2 * k + 3));
    const auto c_o2 = -m_O * k * k / double((2 * k + 1) * (2 * k - 1));
    const auto c_p0 = -m_P * k * (k + 1) / double(2 * (2 * k + 1));

    // "M1" and "O1" (s/X) part:
    const auto cf1 = factor * (c_m1 + c_o1) * (eacp1 + 1.0);
    const auto cg1 = factor * (c_m1 + c_o1) * (eacp1 - 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto s0 = eta_kp1 * (gbk.g0_plus[i] + gbk.gi_plus[i]) +
                      gbk.b0_plus[i] + gbk.bi_plus[i];
      out.f(i) += cf1 * s0 * Fc.g(i);
      out.g(i) += cg1 * s0 * Fc.f(i);
    }

    if (k != 0.0) {
      // "M2" and "O2" (t/Y) part:
      const auto cf2 = factor * (c_m2 + c_o2) * (eac - 1.0);
      const auto cg2 = factor * (c_m2 + c_o2) * (eac + 1.0);
      for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
        const auto t0 = eta_k * (gbk.g0_minus[i] + gbk.gi_minus[i]) -
                        gbk.b0_minus[i] - gbk.bi_minus[i];
        out.f(i) += cf2 * t0 * Fc.g(i);
        out.g(i) += cg2 * t0 * Fc.f(i);
      }

      // "P1" (v/X) part:
      const auto cf3 = factor * c_p0 * (eacp1 + 1.0);
      const auto cg3 = factor * c_p0 * (eacp1 - 1.0);
      for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
        const auto v0 = eta_k * (gbk.g0_minus[i] - gbk.g0_plus[i]) -
                        gbk.b0_minus[i] + gbk.b0_plus[i];
        out.f(i) += cf3 * v0 * Fc.g(i);
        out.g(i) += cg3 * v0 * Fc.f(i);
      }

      // "P2" (w/Y) part:
      const auto cf4 = factor * c_p0 * (eac - 1.0);
      const auto cg4 = factor * c_p0 * (eac + 1.0);
      for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
        const auto w0 = eta_kp1 * (gbk.gi_minus[i] - gbk.gi_plus[i]) +
                        gbk.bi_minus[i] - gbk.bi_plus[i];
        out.f(i) += cf4 * w0 * Fc.g(i);
        out.g(i) += cg4 * w0 * Fc.f(i);
      }
    }
  }

  if (have_n && m_N != 0.0) {
    // "n" part:
    const auto gbk = Breit_gb::single_k_n(Fb, Fd, k);
    const auto factor = k == 0.0 ?
                            0.0 :
                            -m_N * (ka + kc) * (kb + kd) / double(k * (k + 1)) *
                                (sign * m_scale * Dkac * Dkbd);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      out.f(i) += factor * gbk.g[i] * Fc.g(i);
      out.g(i) += factor * gbk.g[i] * Fc.f(i);
    }
  }

  return out;
}

//==============================================================================
DiracSpinor Breit::BPkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                            const DiracSpinor &Fc,
                            const DiracSpinor &Fd) const {

  const auto ka = kappa_v;
  const auto tja = Angular::twoj_k(ka);
  const auto tjb = Fb.twoj();
  const auto tjc = Fc.twoj();
  const auto tjd = Fd.twoj();

  const auto min_twol = std::max(std::abs(tjd - tja), std::abs(tjc - tjb));
  const auto max_twol = std::min(tjd + tja, tjc + tjb);

  DiracSpinor out(0, ka, Fc.grid_sptr());

  for (int twol = min_twol; twol <= max_twol; twol += 2) {
    const auto sj = Angular::sixj_2(tja, tjc, 2 * k, tjb, tjd, twol);
    if (sj == 0.0)
      continue;
    out += sj * Bkv_bcd(twol / 2, kappa_v, Fb, Fd, Fc);
  }
  out *= (2.0 * k + 1.0);
  return out;
}

//==============================================================================
DiracSpinor Breit::dVbrD_Fa(int kappa, int K, const DiracSpinor &Fa,
                            const DiracSpinor &Fb, const DiracSpinor &Xbeta,
                            const DiracSpinor &Ybeta) const {

  const auto twoj = Angular::twoj_k(kappa);
  const auto sign = Angular::neg1pow_2(twoj + Xbeta.twoj() + 2 * K + 1);
  return sign * (BWkv_bcd(K, kappa, Fb, Fa, Xbeta) +
                 BWkv_bcd(K, kappa, Ybeta, Fa, Fb));
}

//==============================================================================
namespace Breit_gb {

void single_k_mop::calculate(const DiracSpinor &Fi, const DiracSpinor &Fj,
                             int k) {

  const auto maxi = std::max(Fi.max_pt(), Fj.max_pt()); // ok?

  // g^{k+1}, g^{k-1}, b^{k+1}, b^{k-1} used in m, o, p
  // Only used if C^k_ij is non-zero
  if (Angular::Ck_kk_SR(k, Fi.kappa(), Fj.kappa())) {
    // assert(k > 0); // ?
    if (k > 0) {
      Coulomb::bk_ab(Fi, Fj, (k - 1), b0_minus, bi_minus, maxi);
      Coulomb::gk_ab(Fi, Fj, (k - 1), g0_minus, gi_minus, maxi);
    }
    Coulomb::bk_ab(Fi, Fj, (k + 1), b0_plus, bi_plus, maxi);
    Coulomb::gk_ab(Fi, Fj, (k + 1), g0_plus, gi_plus, maxi);
  }
}

void single_k_n::calculate(const DiracSpinor &Fi, const DiracSpinor &Fj,
                           int k) {

  const auto maxi = std::max(Fi.max_pt(), Fj.max_pt()); // ok?

  // g^{k} used in n
  // Only used if C^k_-i,j is non-zero
  if (Angular::Ck_kk_SR(k, -Fi.kappa(), Fj.kappa())) {
    Coulomb::gk_ab(Fi, Fj, k, g, gi, maxi);
    using namespace qip::overloads;
    g += gi;
  }
}

} // namespace Breit_gb

} // namespace HF
