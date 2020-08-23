#include "HF/Breit.hpp"
#include "Angular/Angular_369j.hpp"
#include "Coulomb/Coulomb.hpp"

namespace HF {

//******************************************************************************
// Calculates V_br*Fa = \sum_b\sum_k B^k_ba F_b
DiracSpinor Breit::VbrFa(const DiracSpinor &Fa) const {
  DiracSpinor BFa(Fa.n, Fa.k, Fa.rgrid); //
  for (const auto &Fb : *p_core) {
    BkbaFb(&BFa, Fa, Fb);
  }
  return BFa;
}

//******************************************************************************
// Calculates \sum_k B^k_ba F_b
void Breit::BkbaFb(DiracSpinor *BFb, const DiracSpinor &Fa,
                   const DiracSpinor &Fb) const {

  const hidden::Breit_Bk_ba Bkba(Fb, Fa); // every k
  const auto ka = Fa.k;
  const auto kb = Fb.k;

  const auto kmax = (Fb.twoj() + Fa.twoj()) / 2;
  const auto tjap1 = Fa.twoj() + 1;
  for (auto k = 0; k <= kmax; ++k) {
    const auto sk = std::size_t(k);
    const auto Ckba = Angular::Ck_kk(k, Fb.k, Fa.k);

    // M, O, P
    if (Ckba != 0.0) {
      const auto Ckba2 = m_scale * Ckba * Ckba / tjap1;
      const auto [m1, m2] = Mk(k);
      const auto [o1, o2] = Ok(k);
      const auto p1 = Pk(k);

      const auto &b0p = Bkba.bk_0[sk + 1];
      const auto &b0m = Bkba.bk_0[sk - 1]; // XXX
      const auto &bip = Bkba.bk_inf[sk + 1];
      const auto &bim = Bkba.bk_inf[sk - 1]; // XXX
      const auto &g0p = Bkba.gk_0[sk + 1];
      const auto &g0m = Bkba.gk_0[sk - 1]; // XXX
      const auto &gip = Bkba.gk_inf[sk + 1];
      const auto &gim = Bkba.gk_inf[sk - 1]; // XXX
      // XXX -> dangerous: These are never called for k=0, but are INVALID in
      // that case. Even defining the references is undefined behaviour, so
      // this should be fixed. Works though, so..

      const auto ep = eta(k + 1, kb, ka);
      const auto e = eta(k, kb, ka);

      // M1 and O1 term:
      if (m1 != 0.0 || o1 != 0.0) {
        for (auto i = Fb.p0; i < Fb.pinf; ++i) {
          const auto ff =
              Ckba2 * (m1 + o1) * (b0p[i] + bip[i] + ep * g0p[i] + ep * gip[i]);
          BFb->f[i] += ff * (1.0 - ep) * Fb.g[i];
          BFb->g[i] += -ff * (1.0 + ep) * Fb.f[i];
        }
      }

      // M2 and O2 term:
      if (m2 != 0.0 || o2 != 0.0) {
        for (auto i = Fb.p0; i < Fb.pinf; ++i) {
          const auto ff =
              Ckba2 * (m2 + o2) * (-b0m[i] - bim[i] + e * g0m[i] + e * gim[i]);
          BFb->f[i] += -ff * (1.0 + e) * Fb.g[i];
          BFb->g[i] += ff * (1.0 - e) * Fb.f[i];
        }
      }

      // P1 and P2term:
      if (p1 != 0.0) {
        for (auto i = Fb.p0; i < Fb.pinf; ++i) {
          const auto ff =
              Ckba2 * p1 * (-b0m[i] + e * g0m[i] + b0p[i] - e * g0p[i]);
          BFb->f[i] += ff * (1.0 - ep) * Fb.g[i];
          BFb->g[i] += -ff * (1.0 + ep) * Fb.f[i];
        }
        for (auto i = Fb.p0; i < Fb.pinf; ++i) {
          const auto ff =
              // Ckba2 * p1 * (bim[i] + ep * gim[i] - bip[i] + ep * gip[i]);
              Ckba2 * p1 * (bim[i] + ep * gim[i] - bip[i] - ep * gip[i]);
          BFb->f[i] += -ff * (1.0 + e) * Fb.g[i];
          BFb->g[i] += ff * (1.0 - e) * Fb.f[i];
        }
      }

      //
    }

    // N term
    const auto Ckmba = Angular::Ck_kk(k, -Fb.k, Fa.k);
    if (Ckmba != 0.0) {
      const auto Cg = m_scale * Ckmba * Ckmba * Nkba(k, Fb.k, Fa.k) / tjap1;
      const auto &g0 = Bkba.gk_0[sk];
      const auto &gi = Bkba.gk_inf[sk];
      for (auto i = Fb.p0; i < Fb.pinf; ++i) {
        BFb->f[i] += Cg * (g0[i] + gi[i]) * Fb.g[i];
        BFb->g[i] += Cg * (g0[i] + gi[i]) * Fb.f[i];
      }
    }

  } // sum over k
}

//****************************************************************************
// Returns (M+O+P)^k_ba //(for TDHF)
// dV_Br Fa, for given b Beta
// K is multipolarity of RPA operator, k is multipolarity of Breit/Coulomb
DiracSpinor Breit::dVbrFa(int kappa, int K, const DiracSpinor &Fa,
                          const DiracSpinor &Fb, const DiracSpinor &Xbeta,
                          const DiracSpinor &Ybeta) const {

  DiracSpinor dVFa(0, kappa, Fa.rgrid);

  const hidden::Breit_Bk_ba Bkba(Fb, Fa);     // every k
  const hidden::Breit_Bk_ba BkYBa(Ybeta, Fa); // every k
  const auto kn = dVFa.k;
  const auto ka = Fa.k;
  const auto kb = Fb.k;
  const auto kB = Xbeta.k; // kappa(X) = kappa(Y)
  const auto tjn = dVFa.twoj();
  const auto tja = Fa.twoj();
  const auto tjb = Fb.twoj();
  const auto tjB = Xbeta.twoj();

  // (-1)^{jB-ja}
  const auto s_Bb = Angular::evenQ_2(tjB - tja);

  const auto kmax = std::max({tja, tjb, tjB}); // can do better?
  for (auto k = 0; k <= kmax; ++k) {
    const auto s_Kk = Angular::evenQ(K + k);
    const auto sk = std::size_t(k);
    const auto Ckab = Angular::Ck_kk(k, ka, kb);
    const auto Ckab_N = Angular::Ck_kk(k, ka, -kb);
    const auto CknB = Angular::Ck_kk(k, kn, kB);
    const auto sjX = Angular::sixj_2(tja, tjn, 2 * K, tjB, tjb, 2 * k);

    const auto CkaB = Angular::Ck_kk(k, ka, kB);
    const auto Cknb = Angular::Ck_kk(k, kn, kb);
    const auto Cknb_N = Angular::Ck_kk(k, kn, -kb);
    const auto sjY = Angular ::sixj_2(tja, tjn, 2 * K, tjb, tjB, 2 * k);

    const auto cangX = s_Bb * s_Kk * Ckab * CknB * sjX;
    const auto cangY = s_Bb * s_Kk * CkaB * Cknb * sjY;

    const auto cangX_N = m_scale * s_Bb * s_Kk * Ckab_N * CknB * sjX;
    const auto cangY_N = m_scale * s_Bb * s_Kk * CkaB * Cknb_N * sjY;

    // M, O, P, 'X' term
    if (cangX != 0.0) {
      const auto &b0p = Bkba.bk_0[sk + 1];
      const auto &b0m = Bkba.bk_0[sk - 1]; // XXX
      const auto &bip = Bkba.bk_inf[sk + 1];
      const auto &bim = Bkba.bk_inf[sk - 1]; // XXX
      const auto &g0p = Bkba.gk_0[sk + 1];
      const auto &g0m = Bkba.gk_0[sk - 1]; // XXX
      const auto &gip = Bkba.gk_inf[sk + 1];
      const auto &gim = Bkba.gk_inf[sk - 1]; // XXX
      // XXX -> dangerous: These are never called for k=0, but are INVALID in
      // that case. Even defining the references is undefined behaviour, so
      // this should be fixed. Works though, so..
      const auto [m1, m2] = Mk(k);
      const auto [o1, o2] = Ok(k);
      const auto p1 = Pk(k);

      const auto ep = eta(k + 1, kb, ka);
      const auto e = eta(k, kb, ka);

      // M1 and O1 term:
      if (m1 != 0.0 || o1 != 0.0) {
        for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
          const auto ff =
              cangX * (m1 + o1) * (b0p[i] + bip[i] + ep * g0p[i] + ep * gip[i]);
          dVFa.f[i] += ff * (1.0 - ep) * Xbeta.g[i];
          dVFa.g[i] += -ff * (1.0 + ep) * Xbeta.f[i];
        }
      }

      // M2 and O2 term:
      if (m2 != 0.0 || o2 != 0.0) {
        for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
          const auto ff =
              cangX * (m2 + o2) * (-b0m[i] - bim[i] + e * g0m[i] + e * gim[i]);
          dVFa.f[i] += -ff * (1.0 + e) * Xbeta.g[i];
          dVFa.g[i] += ff * (1.0 - e) * Xbeta.f[i];
        }
      }

      // P1 and P2term:
      if (p1 != 0.0) {
        for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
          const auto ff =
              cangX * p1 * (-b0m[i] + e * g0m[i] + b0p[i] - e * g0p[i]);
          dVFa.f[i] += ff * (1.0 - ep) * Xbeta.g[i];
          dVFa.g[i] += -ff * (1.0 + ep) * Xbeta.f[i];
        }
        for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
          const auto ff =
              cangX * p1 * (bim[i] + ep * gim[i] - bip[i] - ep * gip[i]);
          dVFa.f[i] += -ff * (1.0 + e) * Xbeta.g[i];
          dVFa.g[i] += ff * (1.0 - e) * Xbeta.f[i];
        }
      }
    }

    // M, O, P, 'Y' term
    if (cangY != 0.0) {
      const auto &b0p = BkYBa.bk_0[sk + 1];
      const auto &b0m = BkYBa.bk_0[sk - 1]; // XXX
      const auto &bip = BkYBa.bk_inf[sk + 1];
      const auto &bim = BkYBa.bk_inf[sk - 1]; // XXX
      const auto &g0p = BkYBa.gk_0[sk + 1];
      const auto &g0m = BkYBa.gk_0[sk - 1]; // XXX
      const auto &gip = BkYBa.gk_inf[sk + 1];
      const auto &gim = BkYBa.gk_inf[sk - 1]; // XXX
      // XXX -> dangerous: These are never called for k=0, but are INVALID in
      // that case. Even defining the references is undefined behaviour, so
      // this should be fixed. Works though, so..
      const auto [m1, m2] = Mk(k);
      const auto [o1, o2] = Ok(k);
      const auto p1 = Pk(k);

      const auto ep = eta(k + 1, kb, ka);
      const auto e = eta(k, kb, ka);

      // M1 and O1 term:
      if (m1 != 0.0 || o1 != 0.0) {
        for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
          const auto ff =
              cangY * (m1 + o1) * (b0p[i] + bip[i] + ep * g0p[i] + ep * gip[i]);
          dVFa.f[i] += ff * (1.0 - ep) * Xbeta.g[i];
          dVFa.g[i] += -ff * (1.0 + ep) * Xbeta.f[i];
        }
      }

      // M2 and O2 term:
      if (m2 != 0.0 || o2 != 0.0) {
        for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
          const auto ff =
              cangY * (m2 + o2) * (-b0m[i] - bim[i] + e * g0m[i] + e * gim[i]);
          dVFa.f[i] += -ff * (1.0 + e) * Xbeta.g[i];
          dVFa.g[i] += ff * (1.0 - e) * Xbeta.f[i];
        }
      }

      // P1 and P2term:
      if (p1 != 0.0) {
        for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
          const auto ff =
              cangY * p1 * (-b0m[i] + e * g0m[i] + b0p[i] - e * g0p[i]);
          dVFa.f[i] += ff * (1.0 - ep) * Xbeta.g[i];
          dVFa.g[i] += -ff * (1.0 + ep) * Xbeta.f[i];
        }
        for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
          const auto ff =
              cangY * p1 * (bim[i] + ep * gim[i] - bip[i] - ep * gip[i]);
          dVFa.f[i] += -ff * (1.0 + e) * Xbeta.g[i];
          dVFa.g[i] += ff * (1.0 - e) * Xbeta.f[i];
        }
      }
    }

    // N term, 'X'
    if (cangX_N != 0.0) {
      const auto Cg = cangX_N * Nkba(k, Fb.k, Fa.k);
      const auto &g0 = Bkba.gk_0[sk];
      const auto &gi = Bkba.gk_inf[sk];
      for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
        dVFa.f[i] += Cg * (g0[i] + gi[i]) * Xbeta.g[i];
        dVFa.g[i] += Cg * (g0[i] + gi[i]) * Xbeta.f[i];
      }
    }

    // N term, 'Y'
    if (cangY_N != 0.0) {
      const auto Cg = cangY_N * Nkba(k, Fb.k, Fa.k);
      const auto &g0 = BkYBa.gk_0[sk];
      const auto &gi = BkYBa.gk_inf[sk];
      for (auto i = Xbeta.p0; i < Xbeta.pinf; ++i) {
        dVFa.f[i] += Cg * (g0[i] + gi[i]) * Xbeta.g[i];
        dVFa.g[i] += Cg * (g0[i] + gi[i]) * Xbeta.f[i];
      }
    }

  } // sum over k

  return dVFa;
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
    const auto maxi = 0; // std::max(Fa.pinf, Fb.pinf);
    Coulomb::bk_ab(Fb, Fa, int(k), bk_0[k], bk_inf[k], maxi);
    Coulomb::gk_ab(Fb, Fa, int(k), gk_0[k], gk_inf[k], maxi);
  }
}

} // namespace hidden

} // namespace HF
