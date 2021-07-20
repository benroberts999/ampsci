#include "HF/Breit.hpp"
#include "Angular/Angular_369j.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Maths/Grid.hpp"
#include <iostream>
#include <memory>

namespace HF {

//******************************************************************************
// Calculates V_br*Fa = \sum_b\sum_k B^k_ba F_b - For HF potential
DiracSpinor Breit::VbrFa(const DiracSpinor &Fa) const {
  DiracSpinor BFa(Fa.n, Fa.k, Fa.rgrid); //
  for (const auto &Fb : *p_core) {
    BkbaFb(&BFa, Fa, Fb);
  }
  return BFa;
}

//******************************************************************************
// Calculates \sum_k B^k_ba F_b - For HF potential
void Breit::BkbaFb(DiracSpinor *BFb, const DiracSpinor &Fa,
                   const DiracSpinor &Fb) const {
  if (m_scale == 0.0)
    return;

  const hidden::Breit_Bk_ba Bkba(Fb, Fa); // every k
  const auto ka = Fa.k;
  const auto kb = Fb.k;

  const auto kmax = (Fb.twoj() + Fa.twoj()) / 2;
  const auto tjap1 = Fa.twoj() + 1;
  for (auto k = 0; k <= kmax; ++k) {
    const auto Ckba = Angular::Ck_kk(k, kb, ka);

    // M, O, P
    if (Ckba != 0.0) {
      const auto Ckba2 = Ckba * Ckba / tjap1;
      MOPk_ij_Fc(BFb, Ckba2, Bkba, k, kb, ka, Fb);
      //
    }

    // N term
    const auto Ckmba = Angular::Ck_kk(k, -kb, ka);
    if (Ckmba != 0.0) {
      const auto Cg = Ckmba * Ckmba / tjap1;
      Nk_ij_Fc(BFb, Cg, Bkba, k, kb, ka, Fb);
    }

  } // sum over k
}

//******************************************************************************
// Returns (M+O+P)^k_ba //(for TDHF)
// dV_Br Fa, for given b Beta
// K is multipolarity of RPA operator, k is multipolarity of Breit/Coulomb
DiracSpinor Breit::dVbrX_Fa(int kappa, int K, const DiracSpinor &Fa,
                            const DiracSpinor &Fb, const DiracSpinor &Xbeta,
                            const DiracSpinor &Ybeta) const {

  DiracSpinor dVFa(0, kappa, Fa.rgrid);
  // These will be updated in MOPk_ij_Fc if need be
  dVFa.set_min_pt() = Fa.min_pt();
  dVFa.set_max_pt() = Fa.max_pt();
  if (m_scale == 0.0)
    return dVFa;

  // const hidden::Breit_Bk_ba Bkba(Fb, Fa);     // every k
  // const hidden::Breit_Bk_ba BkYBa(Ybeta, Fa); // every k
  // Had to swap these around in order to get agreement.. not sure why?
  // const hidden::Breit_Bk_ba Bkba(Fa, Fb);     // every k
  // const hidden::Breit_Bk_ba BkYBa(Fa, Ybeta); // every k
  // Only create these if at least one angular factor is non-zero
  std::unique_ptr<hidden::Breit_Bk_ba> Bkba{nullptr};
  std::unique_ptr<hidden::Breit_Bk_ba> BkYBa{nullptr};

  const auto kn = dVFa.k;
  const auto ka = Fa.k;
  const auto kb = Fb.k;
  const auto kB = Xbeta.k;
  const auto tjn = dVFa.twoj();
  const auto tja = Fa.twoj();
  const auto tjb = Fb.twoj();
  const auto tjB = Xbeta.twoj();

  // (-1)^{jB-ja}
  const auto s_Bb = Angular::neg1pow_2(tjB - tja);

  const auto kmax = std::max({tja, tjb, tjB}) + K; // ?
  for (auto k = 0; k <= kmax; ++k) {

    const auto s_Kk = Angular::neg1pow(K + k);

    const auto Ckab = Angular::Ck_kk(k, ka, kb);
    const auto Ckab_N = Angular::Ck_kk(k, ka, -kb);
    const auto CknB = Angular::Ck_kk(k, kn, kB);
    const auto sjX = Angular::sixj_2(tja, tjn, 2 * K, tjB, tjb, 2 * k);

    const auto CkaB = Angular::Ck_kk(k, ka, kB);
    const auto Cknb = Angular::Ck_kk(k, kn, kb);
    const auto Cknb_N = Angular::Ck_kk(k, kn, -kb);
    const auto sjY = Angular::sixj_2(tja, tjn, 2 * K, tjb, tjB, 2 * k);

    const auto cangX = s_Bb * s_Kk * Ckab * CknB * sjX;
    const auto cangY = s_Bb * s_Kk * CkaB * Cknb * sjY;
    const auto cangX_N = s_Bb * s_Kk * Ckab_N * CknB * sjX;
    const auto cangY_N = s_Bb * s_Kk * CkaB * Cknb_N * sjY;

    if (!Bkba && (cangX != 0.0 || cangX_N != 0.0)) {
      Bkba = std::make_unique<hidden::Breit_Bk_ba>(Fa, Fb); // every k
    }

    if (!BkYBa && (cangY != 0.0 || cangY_N != 0.0)) {
      // const hidden::Breit_Bk_ba BkYBa(Ybeta, Fa); // every k
      // Had to swap these around in order to get agreement.. not sure why?
      BkYBa = std::make_unique<hidden::Breit_Bk_ba>(Fa, Ybeta); // every k
    }

    if (cangX != 0.0)
      MOPk_ij_Fc(&dVFa, cangX, *Bkba, k, kb, ka, Xbeta);

    if (cangX_N != 0.0)
      Nk_ij_Fc(&dVFa, cangX_N, *Bkba, k, kb, ka, Xbeta);

    if (cangY != 0.0)
      MOPk_ij_Fc(&dVFa, cangY, *BkYBa, k, kB, ka, Fb);

    if (cangY_N != 0.0)
      Nk_ij_Fc(&dVFa, cangY_N, *BkYBa, k, kB, ka, Fb);

  } // sum over k

  return dVFa;
}

//******************************************************************************
DiracSpinor Breit::dVbrD_Fa(int kappa, int K, const DiracSpinor &Fa,
                            const DiracSpinor &Fb, const DiracSpinor &Xbeta,
                            const DiracSpinor &Ybeta) const {

  // XXX Note: This one is hard to check, since it seems to give negligable
  // contribution? Essentially agrees with VD for E1 and E2 (with very tiny
  // assymetry? - so probably minor mistake somewhere!)

  DiracSpinor dVFa(0, kappa, Fa.rgrid);
  // These will be updated in MOPk_ij_Fc if need be
  dVFa.set_min_pt() = Fa.min_pt();
  dVFa.set_max_pt() = Fa.max_pt();
  if (m_scale == 0.0)
    return dVFa;

  // Only need single k here?
  // Also, when X = Y, can do just once!

  const auto kn = dVFa.k;
  const auto ka = Fa.k;
  const auto kb = Fb.k;
  const auto kB = Xbeta.k; // kappa(X) = kappa(Y)

  const auto k = K;

  const auto Ckna = Angular::Ck_kk(k, kn, ka);
  const auto CkBb = Angular::Ck_kk(k, kB, kb);
  const auto CkBb_N = Angular::Ck_kk(k, kB, -kb);

  const auto cang = -Ckna * CkBb / (2 * k + 1);
  const auto cang_N = -Ckna * CkBb_N / (2 * k + 1);

  if (!Angular::zeroQ(cang) || !Angular::zeroQ(cang)) {
    const hidden::Breit_Bk_ba BkXBb(Fb, Xbeta); // every k
    // XXX When X=Y (w=0, real), These are the same(?)!!
    const hidden::Breit_Bk_ba BkYBb(Ybeta, Fb); // every k

    // This "passes" too = check in more detail later!
    // const hidden::Breit_Bk_ba BkYBb(Fb, Ybeta); // every k

    // M, O, P, 'X' term
    if (cang != 0.0) {
      MOPk_ij_Fc(&dVFa, cang, BkXBb, k, kb, kB, Fa); // X
      MOPk_ij_Fc(&dVFa, cang, BkYBb, k, kB, kb, Fa); // Y
    }
    // N term, 'X'
    if (cang_N != 0.0) {
      Nk_ij_Fc(&dVFa, cang_N, BkXBb, k, kb, kB, Fa);
      Nk_ij_Fc(&dVFa, cang_N, BkYBb, k, kB, kb, Fa);
    }
  }

  return dVFa;
}

//******************************************************************************
void Breit::MOPk_ij_Fc(DiracSpinor *BFc, const double Cang,
                       const hidden::Breit_Bk_ba &Bkab, int k, int ki, int kj,
                       const DiracSpinor &Fc) const {
  const auto Cc = m_scale * Cang;
  const auto [m1, m2] = Mk(k);
  const auto [o1, o2] = Ok(k);
  const auto p1 = Pk(k);

  const auto skp = std::size_t(k + 1);
  const auto skm = k != 0 ? std::size_t(k - 1) : 0;
  // nb: when k=0, the 'k-1' vectors are never called
  const auto &b0p = Bkab.bk_0[skp];
  const auto &b0m = Bkab.bk_0[skm];
  const auto &bip = Bkab.bk_inf[skp];
  const auto &bim = Bkab.bk_inf[skm];
  const auto &g0p = Bkab.gk_0[skp];
  const auto &g0m = Bkab.gk_0[skm];
  const auto &gip = Bkab.gk_inf[skp];
  const auto &gim = Bkab.gk_inf[skm];

  const auto ep = eta(k + 1, ki, kj);
  const auto e = eta(k, ki, kj);

  if (Fc.max_pt() > BFc->max_pt())
    BFc->set_max_pt() = Fc.max_pt();
  if (Fc.min_pt() < BFc->min_pt())
    BFc->set_min_pt() = Fc.min_pt();

  // M1 and O1 term:
  if (m1 != 0.0 || o1 != 0.0) {
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto ff =
          Cc * (m1 + o1) * (b0p[i] + bip[i] + ep * g0p[i] + ep * gip[i]);
      BFc->set_f(i) += ff * (1.0 - ep) * Fc.g(i);
      BFc->set_g(i) += -ff * (1.0 + ep) * Fc.f(i);
    }
  }

  // M2 and O2 term:
  if (m2 != 0.0 || o2 != 0.0) {
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto ff =
          Cc * (m2 + o2) * (-b0m[i] - bim[i] + e * g0m[i] + e * gim[i]);
      BFc->set_f(i) += -ff * (1.0 + e) * Fc.g(i);
      BFc->set_g(i) += ff * (1.0 - e) * Fc.f(i);
    }
  }

  // P1 and P2term:
  if (p1 != 0.0) {
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto ff = Cc * p1 * (-b0m[i] + e * g0m[i] + b0p[i] - e * g0p[i]);
      BFc->set_f(i) += ff * (1.0 - ep) * Fc.g(i);
      BFc->set_g(i) += -ff * (1.0 + ep) * Fc.f(i);
    }
    for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
      const auto ff = Cc * p1 * (bim[i] + ep * gim[i] - bip[i] - ep * gip[i]);
      BFc->set_f(i) += -ff * (1.0 + e) * Fc.g(i);
      BFc->set_g(i) += ff * (1.0 - e) * Fc.f(i);
    }
  }
}

//******************************************************************************
void Breit::Nk_ij_Fc(DiracSpinor *BFc, const double Cang,
                     const hidden::Breit_Bk_ba &Bkij, int k, int ki, int kj,
                     const DiracSpinor &Fc) const {

  if (Fc.max_pt() > BFc->max_pt())
    BFc->set_max_pt() = Fc.max_pt();
  if (Fc.min_pt() < BFc->min_pt())
    BFc->set_min_pt() = Fc.min_pt();

  const auto sk = std::size_t(k);
  const auto Cg = m_scale * Cang * Nkba(k, ki, kj);
  const auto &g0 = Bkij.gk_0[sk];
  const auto &gi = Bkij.gk_inf[sk];
  for (auto i = Fc.min_pt(); i < Fc.max_pt(); ++i) {
    BFc->set_f(i) += Cg * (g0[i] + gi[i]) * Fc.g(i);
    BFc->set_g(i) += Cg * (g0[i] + gi[i]) * Fc.f(i);
  }
}

//******************************************************************************
double Breit::eta(int k, int ka, int kb) const {
  return k == 0 ? 0.0 : double(ka - kb) / k;
}

std::pair<double, double> Breit::Mk(int k) const {
  return {-double(k + 1) / (2 * k + 3), -double(k) / (2 * k - 1)};
}

double Breit::Nkba(int k, int kb, int ka) const {
  return k == 0 ? 0.0 : double((kb + ka) * (kb + ka)) / double(k * (k + 1));
}

std::pair<double, double> Breit::Ok(int k) const {
  return {double((k + 1) * (k + 1)) / double((2 * k + 1) * (2 * k + 3)),
          double(k * k) / double((2 * k + 1) * (2 * k - 1))};
}

double Breit::Pk(int k) const {
  return double(k * (k + 1)) / double(2 * (2 * k + 1));
}

//******************************************************************************
namespace hidden {

Breit_Bk_ba::Breit_Bk_ba(const DiracSpinor &Fb, const DiracSpinor &Fa)
    : max_k(std::size_t((Fb.twoj() + Fa.twoj()) / 2) + 1) {

  // size vectors according to kmax :
  bk_0.resize(max_k + 1);
  bk_inf.resize(max_k + 1);
  gk_0.resize(max_k + 1);
  gk_inf.resize(max_k + 1);

  // Fill the b^k functions
  for (auto k = 0ul; k <= max_k; ++k) {
    // a) selection rules (careful)
    // const auto maxi = 0; //
    const auto maxi = std::max(Fa.max_pt(), Fb.max_pt()); // ok?
    // const auto maxi = std::min(Fa.max_pt(), Fb.max_pt()); // XXX OK? No.

    // These normally re-sized in gk_ab()... but not if skip due to SRs
    bk_0[k].resize(Fa.rgrid->num_points());
    bk_inf[k].resize(Fa.rgrid->num_points());
    gk_0[k].resize(Fa.rgrid->num_points());
    gk_inf[k].resize(Fa.rgrid->num_points());

    // bk, in M,N,O, only used if C^(k+/-1)_ba non-zero
    if (Angular::Ck_kk_SR(int(k) - 1, Fb.k, Fa.k) ||
        Angular::Ck_kk_SR(int(k) + 1, Fb.k, Fa.k)) {
      Coulomb::bk_ab(Fb, Fa, int(k), bk_0[k], bk_inf[k], maxi);
    }
    // gk, in M,N,O AND N only used if C^(k+/-1)_ab and C^(k)_(-b,a)
    if (Angular::Ck_kk_SR(int(k) - 1, Fb.k, Fa.k) ||
        Angular::Ck_kk_SR(int(k) + 1, Fb.k, Fa.k) ||
        Angular::Ck_kk_SR(int(k), -Fb.k, Fa.k)) {
      Coulomb::gk_ab(Fb, Fa, int(k), gk_0[k], gk_inf[k], maxi);
    }
  }
}

} // namespace hidden

} // namespace HF
