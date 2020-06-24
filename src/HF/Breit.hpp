#pragma once
#include "Angular/Angular_369j.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <utility>
#include <vector>

namespace HF {

class Breit_Bk_ba { // XXX Make "private"
public:
  Breit_Bk_ba(const DiracSpinor &Fb, const DiracSpinor &Fa)
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

  const std::size_t max_k;
  std::vector<std::vector<double>> bk_0{};
  std::vector<std::vector<double>> bk_inf{};
  std::vector<std::vector<double>> gk_0{};
  std::vector<std::vector<double>> gk_inf{};
};

//******************************************************************************
//! Breit (Hartree-Fock Breit) interaction potential
class Breit {
public:
  //! Contains ptr to core: careful if updating core (e.g., in HF)
  Breit(const std::vector<DiracSpinor> &in_core, double in_scale = 1.0)
      : p_core(&in_core), m_scale(in_scale) {}
  // nb: During HF, must update orbitals each time

  //! Breit(Fa) gives VBr*Fa;
  DiracSpinor operator()(const DiracSpinor &Fa) const { return VbrFa(Fa); }

private:
  const std::vector<DiracSpinor> *const p_core;
  const double m_scale;

  // Calculates V_br*Fa = \sum_b\sum_k B^k_ba F_b
  DiracSpinor VbrFa(const DiracSpinor &Fa) const {
    DiracSpinor BFa(Fa.n, Fa.k, Fa.rgrid); //
    for (const auto &Fb : *p_core) {
      BkbaFb(&BFa, Fa, Fb);
    }
    return BFa;
  }

  // Calculates \sum_k B^k_ba F_b
  void BkbaFb(DiracSpinor *BFb, const DiracSpinor &Fa,
              const DiracSpinor &Fb) const {

    const Breit_Bk_ba Bkba(Fb, Fa); // every k
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
            const auto ff = Ckba2 * (m1 + o1) *
                            (b0p[i] + bip[i] + ep * g0p[i] + ep * gip[i]);
            BFb->f[i] += ff * (1.0 - ep) * Fb.g[i];
            BFb->g[i] += -ff * (1.0 + ep) * Fb.f[i];
          }
        }

        // M2 and O2 term:
        if (m2 != 0.0 || o2 != 0.0) {
          for (auto i = Fb.p0; i < Fb.pinf; ++i) {
            const auto ff = Ckba2 * (m2 + o2) *
                            (-b0m[i] - bim[i] + e * g0m[i] + e * gim[i]);
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

  //----------------------
  double eta(int k, int ka, int kb) const {
    return k == 0 ? 0.0 : double(ka - kb) / k;
  }

  std::pair<double, double> Mk(int k) const {
    return {-double(k + 1) / (2 * k + 3), -double(k) / (2 * k - 1)};
  }

  double Nkba(int k, int kb, int ka) const {
    return k == 0 ? 0.0 : double((kb + ka) * (kb + ka)) / double(k * (k + 1));
  }

  std::pair<double, double> Ok(int k) const {
    return {double((k + 1) * (k + 1)) / double((2 * k + 1) * (2 * k + 3)),
            double(k * k) / double((2 * k + 1) * (2 * k - 1))};
  }

  double Pk(int k) const {
    return double(k * (k + 1)) / double(2 * (2 * k + 1));
  }

public:
  Breit &operator=(const Breit &) = delete;
  Breit(const Breit &) = delete;
  ~Breit() = default;
};

} // namespace HF
