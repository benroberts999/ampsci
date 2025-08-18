#include "HF/Breit.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "Maths/Grid.hpp"
#include "Maths/SphericalBessel.hpp"
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

//-------------------------------------------------------------
DiracSpinor Breit::VbrFa_freqw(const DiracSpinor &Fa,
                               const std::vector<DiracSpinor> &core) const {
  DiracSpinor BFa(Fa.n(), Fa.kappa(), Fa.grid_sptr());
  for (const auto &Fb : core) {
    const auto kmin = std::abs(Fb.twoj() - Fa.twoj()) / 2;
    const auto kmax = (Fb.twoj() + Fa.twoj()) / 2;

    const auto w = PhysConst::alpha * std::abs(Fa.en() - Fb.en());

    // if w is too small then the frequency-dependent integrals diverge, so use the frequency-independent equations instead
    if (w < PhysConst::alpha2) {
      for (int k = kmin; k <= kmax; ++k) {
        const auto s = Angular::neg1pow(k);
        if (s == 1) {
          BFa += Bkv_bcd(k, Fa.kappa(), Fb, Fb, Fa);
        } else {
          BFa -= Bkv_bcd(k, Fa.kappa(), Fb, Fb, Fa);
        }
      }
    } else {
      for (int k = kmin; k <= kmax; ++k) {
        const auto s = Angular::neg1pow(k);
        if (s == 1) {
          BFa += Bkv_bcd_freqw(k, Fa.kappa(), Fb, Fb, Fa, w);
        } else {
          BFa -= Bkv_bcd_freqw(k, Fa.kappa(), Fb, Fb, Fa, w);
        }
      }
    }
  }
  BFa *= (-1.0 / Fa.twojp1());
  return BFa;
}

DiracSpinor
Breit::VbrFa_freqw_hermitian(const DiracSpinor &Fa,
                             const std::vector<DiracSpinor> &core,
                             const std::vector<DiracSpinor> &basis) const {

  // initialise spinor to store rsult
  DiracSpinor BFa(Fa.n(), Fa.kappa(), Fa.grid_sptr());

  // sum over core states and complete set of basis states
  for (const auto &Fb : core) {
    const auto kmin = std::abs(Fb.twoj() - Fa.twoj()) / 2;
    const auto kmax = (Fb.twoj() + Fa.twoj()) / 2;
    for (const auto &Fi : basis) {
      for (int k = kmin; k <= kmax; ++k) {
        const auto s = Angular::neg1pow(k);
        if (s == 1) {
          BFa += Bk_abcd(k, Fi, Fb, Fa, Fb) * Fi;
          BFa += Bk_abcd_eac_freqw(k, Fi, Fb, Fa, Fb) * Fi;
          BFa -= Bk_abcd_ebd_freqw(k, Fi, Fb, Fb, Fa) * Fi;
          BFa -= Bk_abcd_eac_freqw(k, Fi, Fb, Fb, Fa) * Fi;
        } else {
          BFa -= Bk_abcd(k, Fi, Fb, Fa, Fb) * Fi;
          BFa -= Bk_abcd_eac_freqw(k, Fi, Fb, Fa, Fb) * Fi;
          BFa += Bk_abcd_ebd_freqw(k, Fi, Fb, Fb, Fa) * Fi;
          BFa += Bk_abcd_eac_freqw(k, Fi, Fb, Fb, Fa) * Fi;
        }
      }
    }
  }
  BFa *= 0.5;
  // multiply by angular factor out the front but is this needed here?
  BFa *= (-1.0 / Fa.twojp1());
  return BFa;
}

//==============================================================================
// can i do something like this for frequency-dependent Breit? The m_gb and m_gb_N vectors will need another index for the energy of the transition - but we could have any arbitrary orbital that we are acting VBr(w) on so this won't work i think
void Breit::fill_gb(const std::vector<DiracSpinor> &basis, int t_max_k) {

  const auto max_k = std::min(DiracSpinor::max_tj(basis), t_max_k);
  m_gb.resize(std::size_t(max_k + 1));
  m_gb_N.resize(std::size_t(max_k + 1));

  for (int k = 1; k <= max_k; ++k) {
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        m_gb.at(std::size_t(k)).add(a, b, Breit_gb::single_k_mop());
        m_gb_N.at(std::size_t(k)).add(a, b, Breit_gb::single_k_n());
      }
    }
  }

#pragma omp parallel for collapse(3)
  for (int k = 1; k <= max_k; ++k) {
    for (std::size_t ia = 0; ia < basis.size(); ++ia) {
      for (std::size_t ib = 0; ib < basis.size(); ++ib) {
        const auto &a = basis.at(ia);
        const auto &b = basis.at(ib);
        m_gb.at(std::size_t(k)).update(a, b, Breit_gb::single_k_mop(a, b, k));
        m_gb_N.at(std::size_t(k)).update(a, b, Breit_gb::single_k_n(a, b, k));
      }
    }
  }
}

//==============================================================================
DiracSpinor Breit::Bkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                           const DiracSpinor &Fc, const DiracSpinor &Fd) const {

  DiracSpinor out(0, kappa_v, Fc.grid_sptr());
  if (k == 0) {
    out.min_pt() = 0;
    out.max_pt() = 0;
    return out;
  }
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

  const auto have_mop = (Ckac != 0.0) && (Ckbd != 0.0);
  const auto have_n =
      (Dkac != 0.0) && (Dkbd != 0.0) && (ka + kc != 0) && (kb + kd != 0);

  // nb: never have both MOP _and_ N ! different parity rule for each k!
  assert(!(have_mop && have_n));

  if (have_mop) {
    // Angular factors
    const auto d_ac = ka - kc;
    const auto d_bd = kb - kd;
    const auto d_bd_k = d_bd / double(k);
    const auto d_bd_kp1 = d_bd / double(k + 1);
    const auto d_ac_kp1 = d_ac / double(k + 1);
    const auto d_ac_k = d_ac / double(k);

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
    const auto cf1 = factor * (c_m1 + c_o1) * (d_ac_kp1 + 1.0);
    const auto cg1 = factor * (c_m1 + c_o1) * (d_ac_kp1 - 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto s0 = d_bd_kp1 * (gbk.g0_plus[i] + gbk.gi_plus[i]) +
                      gbk.b0_plus[i] + gbk.bi_plus[i];
      out.f(i) += cf1 * s0 * Fc.g(i);
      out.g(i) += cg1 * s0 * Fc.f(i);
    }

    // "M2" and "O2" (t/Y) part:
    const auto cf2 = factor * (c_m2 + c_o2) * (d_ac_k - 1.0);
    const auto cg2 = factor * (c_m2 + c_o2) * (d_ac_k + 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto t0 = d_bd_k * (gbk.g0_minus[i] + gbk.gi_minus[i]) -
                      gbk.b0_minus[i] - gbk.bi_minus[i];
      out.f(i) += cf2 * t0 * Fc.g(i);
      out.g(i) += cg2 * t0 * Fc.f(i);
    }

    // "P1" (v/X) part:
    const auto cf3 = factor * c_p0 * (d_ac_kp1 + 1.0);
    const auto cg3 = factor * c_p0 * (d_ac_kp1 - 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto v0 = d_bd_k * (gbk.g0_minus[i] - gbk.g0_plus[i]) -
                      gbk.b0_minus[i] + gbk.b0_plus[i];
      //out.f(i) += cf3 * v0 * Fc.g(i);
      //out.g(i) += cg3 * v0 * Fc.f(i);
    }

    // "P2" (w/Y) part:
    const auto cf4 = factor * c_p0 * (d_ac_k - 1.0);
    const auto cg4 = factor * c_p0 * (d_ac_k + 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto w0 = d_bd_kp1 * (gbk.gi_minus[i] - gbk.gi_plus[i]) +
                      gbk.bi_minus[i] - gbk.bi_plus[i];
      //out.f(i) += cf4 * w0 * Fc.g(i);
      //out.g(i) += cg4 * w0 * Fc.f(i);
    }
  }

  if (have_n && m_N != 0.0) {
    // "n" part:
    const auto gbk = Breit_gb::single_k_n(Fb, Fd, k);
    const auto factor = -m_N * (ka + kc) * (kb + kd) / double(k * (k + 1)) *
                        (sign * m_scale * Dkac * Dkbd);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      out.f(i) += factor * gbk.g[i] * Fc.g(i);
      out.g(i) += factor * gbk.g[i] * Fc.f(i);
    }
  }

  return out;
}

// Calculates (m^k(w) + n^k(w) + o^k(w) + p^k(w))Fi(r) -- frequency-dependent Breit
DiracSpinor Breit::Bkv_bcd_freqw(int k, int kappa_v, const DiracSpinor &Fb,
                                 const DiracSpinor &Fc, const DiracSpinor &Fd,
                                 const double w) const {

  // if energy is zero then just use the frequency-independent Breit HF operator to avoid the Bessel functions being called at zero
  if (w == 0) {
    return Bkv_bcd(k, kappa_v, Fb, Fc, Fd);
  }

  DiracSpinor out(0, kappa_v, Fc.grid_sptr());
  if (k == 0) {
    out.min_pt() = 0;
    out.max_pt() = 0;
    return out;
  }
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

  const auto have_mop = (Ckac != 0.0) && (Ckbd != 0.0);
  const auto have_n =
      (Dkac != 0.0) && (Dkbd != 0.0) && (ka + kc != 0) && (kb + kd != 0);

  // nb: never have both MOP _and_ N ! different parity rule for each k!
  assert(!(have_mop && have_n));

  if (have_mop) {
    // Angular factors
    const auto d_ac = ka - kc;
    const auto d_bd = kb - kd;
    const auto d_bd_k = d_bd / double(k);
    const auto d_bd_kp1 = d_bd / double(k + 1);
    const auto d_ac_kp1 = d_ac / double(k + 1);
    const auto d_ac_k = d_ac / double(k);

    // Calculate the frequency-dependent Breit radial screening integrals
    const auto ghkw = Breit_gh_freqdep::single_k_mop_freq(Fb, Fd, k, w);

    // coeficients (including scaling factors)
    const auto factor = sign * m_scale * Ckac * Ckbd;
    const auto c_m1 = m_M * (k + 1) / double(2 * k + 3);
    const auto c_m2 = m_M * k / double(2 * k - 1);
    const auto c_o1 =
        -m_O * (k + 1) * (k + 1) / double((2 * k + 1) * (2 * k + 3));
    const auto c_o2 = -m_O * k * k / double((2 * k + 1) * (2 * k - 1));
    const auto c_p0 = -m_P * k * (k + 1) / double(2 * (2 * k + 1));

    // "M1" and "O1" (s/X) part:
    const auto cf1 = factor * (c_m1 + c_o1) * (d_ac_kp1 + 1.0);
    const auto cg1 = factor * (c_m1 + c_o1) * (d_ac_kp1 - 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto s0 =
          d_bd_kp1 * (ghkw.g0_plus_freqw[i] + ghkw.gi_plus_freqw[i]) +
          ghkw.h0_plus_freqw[i] + ghkw.hi_plus_freqw[i];
      out.f(i) += cf1 * s0 * Fc.g(i);
      out.g(i) += cg1 * s0 * Fc.f(i);
    }

    // "M2" and "O2" (t/Y) part:
    const auto cf2 = factor * (c_m2 + c_o2) * (d_ac_k - 1.0);
    const auto cg2 = factor * (c_m2 + c_o2) * (d_ac_k + 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto t0 =
          d_bd_k * (ghkw.g0_minus_freqw[i] + ghkw.gi_minus_freqw[i]) -
          ghkw.h0_minus_freqw[i] - ghkw.hi_minus_freqw[i];
      out.f(i) += cf2 * t0 * Fc.g(i);
      out.g(i) += cg2 * t0 * Fc.f(i);
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS TERM IS A PROBLEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // "P1" (v/X) part:
    const auto cf3 = factor * c_p0 * (d_ac_kp1 + 1.0);
    const auto cg3 = factor * c_p0 * (d_ac_kp1 - 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto v1pv4 = ghkw.v1_freqw[i] + ghkw.v4_freqw[i];
      //out.f(i) += cf3 * v1pv4 * Fc.g(i);
      //out.g(i) += cg3 * v1pv4 * Fc.f(i);
    }

    // "P2" (w/Y) part:
    const auto cf4 = factor * c_p0 * (d_ac_k + 1.0);
    const auto cg4 = factor * c_p0 * (d_ac_k - 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto v2pv3 = ghkw.v2_freqw[i] + ghkw.v3_freqw[i];
      //out.f(i) += cf4 * v2pv3 * Fc.g(i);
      //out.g(i) += cg4 * v2pv3 * Fc.f(i);
    }
  }

  if (have_n && m_N != 0.0) {
    // "n" part:
    const auto ghkw = Breit_gh_freqdep::single_k_n_freq(Fb, Fd, k, w);
    const auto factor = -m_N * (ka + kc) * (kb + kd) / double(k * (k + 1)) *
                        (sign * m_scale * Dkac * Dkbd);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      out.f(i) += factor * ghkw.g[i] * Fc.g(i);
      out.g(i) += factor * ghkw.g[i] * Fc.f(i);
    }
  }

  return out;
}

//==============================================================================
double Breit::Bk_abcd_2(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                        const DiracSpinor &Fc, const DiracSpinor &Fd) const {

  assert(!m_gb.empty());

  double out = 0.0;

  if (k == 0) {
    return out;
  }

  const auto min = std::max(Fa.min_pt(), Fc.min_pt());
  const auto max = std::min(Fa.max_pt(), Fc.max_pt());
  const auto &drdu = Fa.grid().drdu();

  const auto ka = Fa.kappa();
  const auto kb = Fb.kappa();
  const auto kc = Fc.kappa();
  const auto kd = Fd.kappa();

  const auto Ckac = Angular::Ck_kk(k, ka, kc);
  const auto Ckbd = Angular::Ck_kk(k, kb, kd);
  const auto Dkac = Angular::Ck_kk(k, -ka, kc);
  const auto Dkbd = Angular::Ck_kk(k, -kb, kd);
  const auto tja = Angular::twoj_k(ka);
  const auto sign = Angular::neg1pow_2(2 * k + tja - Fb.twoj());

  const auto have_mop = (Ckac != 0.0) && (Ckbd != 0.0);
  const auto have_n =
      (Dkac != 0.0) && (Dkbd != 0.0) && (ka + kc != 0) && (kb + kd != 0);

  // nb: never have both MOP _and_ N ! different parity rule for each k!
  assert(!(have_mop && have_n));

  if (have_mop && m_gb.size() > std::size_t(k) &&
      m_gb.at(std::size_t(k)).contains(Fb, Fd)) {
    // Angular factors
    const auto d_ac = ka - kc;
    const auto d_bd = kb - kd;
    const auto d_bd_k = d_bd / double(k);
    const auto d_bd_kp1 = d_bd / double(k + 1);
    const auto d_ac_kp1 = d_ac / double(k + 1);
    const auto d_ac_k = d_ac / double(k);

    // Calculate the Breit radial screening integrals
    // const auto gbk = Breit_gb::single_k_mop(Fb, Fd, k); // this: faster!
    const auto &gbk = m_gb.at(std::size_t(k)).getv(Fb, Fd); // this: faster!

    // coeficients (including scaling factors)
    const auto factor = sign * m_scale * Ckac * Ckbd;
    const auto c_m1 = m_M * (k + 1) / double(2 * k + 3);
    const auto c_m2 = m_M * k / double(2 * k - 1);
    const auto c_o1 =
        -m_O * (k + 1) * (k + 1) / double((2 * k + 1) * (2 * k + 3));
    const auto c_o2 = -m_O * k * k / double((2 * k + 1) * (2 * k - 1));
    const auto c_p0 = -m_P * k * (k + 1) / double(2 * (2 * k + 1));

    // "M1" and "O1" (s/X) part:
    const auto cf1 = factor * (c_m1 + c_o1) * (d_ac_kp1 + 1.0);
    const auto cg1 = factor * (c_m1 + c_o1) * (d_ac_kp1 - 1.0);
    for (auto i = min; i < max; ++i) {
      const auto s0 = d_bd_kp1 * (gbk.g0_plus[i] + gbk.gi_plus[i]) +
                      gbk.b0_plus[i] + gbk.bi_plus[i];
      out += Fa.f(i) * cf1 * s0 * Fc.g(i) * drdu[i];
      out += Fa.g(i) * cg1 * s0 * Fc.f(i) * drdu[i];
    }

    // "M2" and "O2" (t/Y) part:
    const auto cf2 = factor * (c_m2 + c_o2) * (d_ac_k - 1.0);
    const auto cg2 = factor * (c_m2 + c_o2) * (d_ac_k + 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto t0 = d_bd_k * (gbk.g0_minus[i] + gbk.gi_minus[i]) -
                      gbk.b0_minus[i] - gbk.bi_minus[i];
      out += Fa.f(i) * cf2 * t0 * Fc.g(i) * drdu[i];
      out += Fa.g(i) * cg2 * t0 * Fc.f(i) * drdu[i];
    }

    // "P1" (v/X) part:
    const auto cf3 = factor * c_p0 * (d_ac_kp1 + 1.0);
    const auto cg3 = factor * c_p0 * (d_ac_kp1 - 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto v0 = d_bd_k * (gbk.g0_minus[i] - gbk.g0_plus[i]) -
                      gbk.b0_minus[i] + gbk.b0_plus[i];
      out += Fa.f(i) * cf3 * v0 * Fc.g(i) * drdu[i];
      out += Fa.g(i) * cg3 * v0 * Fc.f(i) * drdu[i];
    }

    // "P2" (w/Y) part:
    const auto cf4 = factor * c_p0 * (d_ac_k - 1.0);
    const auto cg4 = factor * c_p0 * (d_ac_k + 1.0);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto w0 = d_bd_kp1 * (gbk.gi_minus[i] - gbk.gi_plus[i]) +
                      gbk.bi_minus[i] - gbk.bi_plus[i];
      out += Fa.f(i) * cf4 * w0 * Fc.g(i) * drdu[i];
      out += Fa.g(i) * cg4 * w0 * Fc.f(i) * drdu[i];
    }
  }

  if (have_n && m_N != 0.0 && m_gb_N.size() > std::size_t(k) &&
      m_gb_N.at(std::size_t(k)).contains(Fb, Fd)) {
    // "n" part:
    // const auto gbk = Breit_gb::single_k_n(Fb, Fd, k);
    const auto &gbk = m_gb_N.at(std::size_t(k)).getv(Fb, Fd);
    const auto factor = -m_N * (ka + kc) * (kb + kd) / double(k * (k + 1)) *
                        (sign * m_scale * Dkac * Dkbd);
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      out += Fa.f(i) * factor * gbk.g[i] * Fc.g(i) * drdu[i];
      out += Fa.g(i) * factor * gbk.g[i] * Fc.f(i) * drdu[i];
    }
  }

  return out * Fa.grid().du();
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

  DiracSpinor out(0, ka, Fd.grid_sptr());
  out.min_pt() = Fd.min_pt();
  out.max_pt() = Fd.max_pt();

  for (int twol = min_twol; twol <= max_twol; twol += 2) {
    const auto sj = Angular::sixj_2(tjc, tja, 2 * k, tjd, tjb, twol);
    if (sj == 0.0)
      continue;
    out += sj * Bkv_bcd(twol / 2, kappa_v, Fb, Fd, Fc);
  }
  out *= (2.0 * k + 1.0);
  return out;
}

//==============================================================================
double Breit::BPk_abcd_2(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                         const DiracSpinor &Fc, const DiracSpinor &Fd) const {

  const auto tja = Fa.twoj();
  const auto tjb = Fb.twoj();
  const auto tjc = Fc.twoj();
  const auto tjd = Fd.twoj();

  const auto min_twol = std::max(std::abs(tjd - tja), std::abs(tjc - tjb));
  const auto max_twol = std::min(tjd + tja, tjc + tjb);

  double out = 0.0;
  for (int twol = min_twol; twol <= max_twol; twol += 2) {
    const auto sj = Angular::sixj_2(tjc, tja, 2 * k, tjd, tjb, twol);
    if (sj == 0.0)
      continue;
    out += sj * Bk_abcd(twol / 2, Fa, Fb, Fd, Fc);
  }
  out *= (2.0 * k + 1.0);
  return out;
}

//==============================================================================
DiracSpinor Breit::BWkv_bcd(int k, int kappa_v, const DiracSpinor &Fb,
                            const DiracSpinor &Fc,
                            const DiracSpinor &Fd) const {
  return Bkv_bcd(k, kappa_v, Fb, Fc, Fd) + BPkv_bcd(k, kappa_v, Fb, Fc, Fd);
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

double Breit::BWk_abcd_2(int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                         const DiracSpinor &Fc, const DiracSpinor &Fd) const {
  return Bk_abcd_2(k, Fa, Fb, Fc, Fd) + BPk_abcd_2(k, Fa, Fb, Fc, Fd);
}

// calculates two-particle Breit matrix element evaluated at the energy difference of electrons a and c
double Breit::Bk_abcd_eac_freqw(int k, const DiracSpinor &Fa,
                                const DiracSpinor &Fb, const DiracSpinor &Fc,
                                const DiracSpinor &Fd) const {
  return Fa * Bkv_bcd_freqw(k, Fa.kappa(), Fb, Fc, Fd,
                            PhysConst::alpha * abs(Fa.en() - Fc.en()));

  // for testing
  //return Fa * Bkv_bcd_freqw(k, Fa.kappa(), Fb, Fc, Fd, PhysConst::alpha2);
}

// calculates two-particle Breit matrix element evaluated at the energy difference of electrons b and d
double Breit::Bk_abcd_ebd_freqw(int k, const DiracSpinor &Fa,
                                const DiracSpinor &Fb, const DiracSpinor &Fc,
                                const DiracSpinor &Fd) const {
  return Fa * Bkv_bcd_freqw(k, Fa.kappa(), Fb, Fc, Fd,
                            PhysConst::alpha * abs(Fb.en() - Fd.en()));

  // for testing
  //return Fa * Bkv_bcd_freqw(k, Fa.kappa(), Fb, Fc, Fd, 1000.0);
}

//==============================================================================
DiracSpinor Breit::dV_Br(int kappa, int K, const DiracSpinor &Fa,
                         const DiracSpinor &Fb, const DiracSpinor &Xbeta,
                         const DiracSpinor &Ybeta) const {

  const auto twoj = Angular::twoj_k(kappa);
  const auto sign = Angular::neg1pow_2(twoj + Xbeta.twoj() + 2 * K + 2);
  const double tkp1 = double(2 * K + 1);
  return (sign / tkp1) * (BWkv_bcd(K, kappa, Fb, Fa, Xbeta) +
                          BWkv_bcd(K, kappa, Ybeta, Fa, Fb));
}

//==============================================================================
double Breit::de2_HF(const DiracSpinor &v,
                     const std::vector<DiracSpinor> &holes,
                     const std::vector<DiracSpinor> &excited) const {
  double deHF = 0.0;
#pragma omp parallel for reduction(+ : deHF)
  for (std::size_t im = 0; im < excited.size(); ++im) {
    const auto &m = excited[im];
    for (const auto &a : holes) {
      if (a.kappa() != m.kappa())
        continue;
      const auto s1 = Angular::neg1pow_2(v.twoj() - m.twoj());
      // .... Would be faster to use one of the tables ....
      deHF -= s1 * Coulomb::Wk_abcd(0, v, m, v, a) * (a * VbrFa(m, holes)) /
              (m.en() - a.en());
      // nb: W symmetic since v=w : W_wavm = W_vmwa
      const auto s2 = Angular::neg1pow_2(v.twoj() - a.twoj());
      deHF -= s2 * Coulomb::Wk_abcd(0, v, a, v, m) * (m * VbrFa(a, holes)) /
              (m.en() - a.en());
    }
  }
  return deHF;
}

double Breit::de2(const DiracSpinor &v, const std::vector<DiracSpinor> &holes,
                  const std::vector<DiracSpinor> &excited) const {
  double de = 0.0;
#pragma omp parallel for reduction(+ : de)
  for (std::size_t in = 0; in < excited.size(); ++in) {
    const auto &n = excited[in];
    for (const auto &a : holes) {
      for (const auto &b : holes) {
        const auto [kmin, kmax] = Coulomb::k_minmax_W(v, n, a, b);
        for (int k = kmin; k <= kmax; k++) {
          const auto denom = v.en() + n.en() - a.en() - b.en();
          const auto f = 2.0 / ((2 * k + 1) * v.twojp1() * denom);
          de += f * Coulomb::Wk_abcd(k, v, n, a, b) * Bk_abcd(k, v, n, a, b);
        }
      }
      for (const auto &m : excited) {
        const auto [kmin, kmax] = Coulomb::k_minmax_W(m, n, v, a);
        for (int k = kmin; k <= kmax; k++) {
          const auto denom = m.en() + n.en() - v.en() - a.en();
          const auto f = 2.0 / ((2 * k + 1) * v.twojp1() * denom);
          de += -f * Coulomb::Wk_abcd(k, m, n, v, a) * Bk_abcd(k, m, n, v, a);
        }
      }
    }
  }
  return de;
}

//==============================================================================
namespace Breit_gb {

void single_k_mop::calculate(const DiracSpinor &Fi, const DiracSpinor &Fj,
                             int k) {

  if (!Angular::Ck_kk_SR(k, Fi.kappa(), Fj.kappa()))
    return;

  const auto maxi =
      std::min({Fi.max_pt(), Fj.max_pt(), Fi.grid().num_points()}); // ok?

  // g^{k+1}, g^{k-1}, b^{k+1}, b^{k-1} used in m, o, p
  assert(k != 0);

#pragma omp parallel sections num_threads(4)
  {
#pragma omp section
    Coulomb::bk_ab((k - 1), Fi, Fj, b0_minus, bi_minus, maxi);
#pragma omp section
    Coulomb::gk_ab((k - 1), Fi, Fj, g0_minus, gi_minus, maxi);
#pragma omp section
    Coulomb::bk_ab((k + 1), Fi, Fj, b0_plus, bi_plus, maxi);
#pragma omp section
    Coulomb::gk_ab((k + 1), Fi, Fj, g0_plus, gi_plus, maxi);
  }
}

void single_k_n::calculate(const DiracSpinor &Fi, const DiracSpinor &Fj,
                           int k) {
  if (!Angular::Ck_kk_SR(k, -Fi.kappa(), Fj.kappa()))
    return;
  const auto maxi =
      std::min({Fi.max_pt(), Fj.max_pt(), Fi.grid().num_points()}); // ok?
  assert(k != 0);
  Coulomb::gk_ab(k, Fi, Fj, g, gi, maxi);
  using namespace qip::overloads;
  g += gi;
}

} // namespace Breit_gb

namespace Breit_gh_freqdep {
// generates frequency-dependent Breit screening functions for calculating m, o, p
void single_k_mop_freq::calculate(const DiracSpinor &Fi, const DiracSpinor &Fj,
                                  int k, const double w) {

  if (!Angular::Ck_kk_SR(k, Fi.kappa(), Fj.kappa()))
    return;

  const auto maxi =
      std::min({Fi.max_pt(), Fj.max_pt(), Fi.grid().num_points()}); // ok?

  // g^{k+1}(w), g^{k-1}(w), h^{k+1}(w), h^{k-1}(w) used in m, o, p
  assert(k != 0);

  const auto r = Fi.grid().r();

// calculates
#pragma omp parallel sections num_threads(4)
  {
#pragma omp section
    Coulomb::gk_ab_freqw(k - 1, Fi, Fj, g0_minus_freqw, gi_minus_freqw, maxi,
                         w);
#pragma omp section
    Coulomb::hk_ab_freqw(k - 1, Fi, Fj, h0_minus_freqw, hi_minus_freqw, maxi,
                         w);
#pragma omp section
    Coulomb::gk_ab_freqw(k + 1, Fi, Fj, g0_plus_freqw, gi_plus_freqw, maxi, w);
#pragma omp section
    Coulomb::hk_ab_freqw(k + 1, Fi, Fj, h0_plus_freqw, hi_plus_freqw, maxi, w);
#pragma omp section
    Coulomb::vk_ab_freqw(k, Fi, Fj, Fi.grid(), v1_freqw, v2_freqw, v3_freqw,
                         v4_freqw, maxi, w);
  }
}

// constructs full g^k - only used for calculating u^k_{abcd} and then from there n^k_{abcd}
void single_k_n_freq::calculate(const DiracSpinor &Fi, const DiracSpinor &Fj,
                                int k, const double w) {
  if (!Angular::Ck_kk_SR(k, -Fi.kappa(), Fj.kappa()))
    return;
  const auto maxi =
      std::min({Fi.max_pt(), Fj.max_pt(), Fi.grid().num_points()}); // ok?
  assert(k != 0);
  Coulomb::gk_ab_freqw(k, Fi, Fj, g, gi, maxi, w);
  using namespace qip::overloads;
  g += gi;
}

} // namespace Breit_gh_freqdep

} // namespace HF
