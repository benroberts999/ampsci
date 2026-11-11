#include "DCP.hpp"
#include "Angular/SixJTable.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cmath>
#include <cstddef>

namespace Amplitudes {

//==============================================================================

double dcp_C1(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
              double omega_s, int K, int kt, int ks,
              const Coulomb::meTable<double> &T,
              const Coulomb::meTable<double> &S,
              const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited) {

  const auto s_Kts = Angular::neg1pow(K + kt + ks);

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak) collapse(2)
  for (std::size_t im = 0; im < excited.size(); ++im) {
    for (std::size_t in = 0; in < excited.size(); ++in) {
      const auto &m = excited[im];
      const auto &n = excited[in];

      const auto s_wm = Angular::neg1pow_2(w.twoj() + m.twoj());
      const auto Wmwvn = Coulomb::Wk_abcd(K, w, m, v, n);
      if (Wmwvn == 0.0)
        continue;

      for (const auto &a : core) {

        // C1_ts:
        const auto sj_ts = Angular::SixJ(kt, ks, K, m, n, a);
        if (sj_ts != 0.0) {
          const auto Tna = T.getv(n, a);
          const auto Sam = S.getv(a, m);
          const auto de_ts =
            (a.en() - n.en() + omega_t) * (a.en() - m.en() - omega_s);
          Ak += Tna * Sam * Wmwvn / de_ts * s_wm * sj_ts;
        }

        // C1_st:
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, m, n, a);
        if (sj_st != 0.0) {
          const auto Sna = S.getv(n, a);
          const auto Tam = T.getv(a, m);
          const auto de_st =
            (a.en() - n.en() + omega_s) * (a.en() - m.en() - omega_t);
          Ak += Sna * Tam * Wmwvn / de_st * s_wm * s_Kts * sj_st;
        }
      }
    }
  }
  return Ak / std::sqrt(2 * K + 1);
}

//==============================================================================

double dcp_C2(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
              double omega_s, int K, int kt, int ks,
              const Coulomb::meTable<double> &T,
              const Coulomb::meTable<double> &S,
              const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited) {

  const auto s_Kts = Angular::neg1pow(K + kt + ks);

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak) collapse(2)
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    for (std::size_t ib = 0; ib < core.size(); ++ib) {
      const auto &a = core[ia];
      const auto &b = core[ib];

      const auto s_wa = Angular::neg1pow_2(w.twoj() - a.twoj());
      const auto Wwavb = Coulomb::Wk_abcd(K, w, a, v, b);
      if (Wwavb == 0.0)
        continue;

      for (const auto &n : excited) {

        // C2_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, b, a, n);
        if (sj_ts != 0.0) {
          const auto Tna = T.getv(n, a);
          const auto Sbn = S.getv(b, n);
          const auto de_ts =
            (a.en() - n.en() + omega_t) * (b.en() - n.en() - omega_s);
          Ak += s_wa * s_Kts * sj_ts * Tna * Sbn * Wwavb / de_ts;
        }

        // C2_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, b, a, n);
        if (sj_st != 0.0) {
          const auto Sna = S.getv(n, a);
          const auto Tbn = T.getv(b, n);
          const auto de_st =
            (a.en() - n.en() + omega_s) * (b.en() - n.en() - omega_t);
          Ak += s_wa * sj_st * Sna * Tbn * Wwavb / de_st;
        }
      }
    }
  }
  return Ak / std::sqrt(2 * K + 1);
}

//==============================================================================

double dcp_R12(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
               double omega_s, int K, int kt, int ks,
               const Coulomb::meTable<double> &T,
               const Coulomb::meTable<double> &S,
               const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited) {

  // Both R1 and R2 share a Wk integral; faster to combine

  const auto s_Kts = Angular::neg1pow(K + kt + ks);
  const auto omega_ts = omega_t + omega_s;

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak) collapse(2)
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    for (std::size_t in = 0; in < excited.size(); ++in) {
      const auto &a = core[ia];
      const auto &n = excited[in];

      const auto Wwavn = Coulomb::Wk_abcd(K, w, a, v, n);
      if (Wwavn == 0.0)
        continue;

      // R1:
      const auto s_wa = Angular::neg1pow_2(w.twoj() + a.twoj());
      for (const auto &m : excited) {

        //R1_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, n, a, m);
        if (sj_ts != 0.0) {
          const auto Tma = T.getv(m, a);
          const auto Snm = S.getv(n, m);
          const auto de_ts =
            (a.en() - m.en() + omega_t) * (a.en() - n.en() + omega_ts);
          Ak += s_wa * s_Kts * sj_ts * Tma * Snm * Wwavn / de_ts;
        }

        //R1_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, n, a, m);
        if (sj_st != 0.0) {
          const auto Sma = S.getv(m, a);
          const auto Tnm = T.getv(n, m);
          const auto de_st =
            (a.en() - m.en() + omega_s) * (a.en() - n.en() + omega_ts);
          Ak += s_wa * sj_st * Sma * Tnm * Wwavn / de_st;
        }
      }

      // R2 (note -ve sign cf R1):
      for (const auto &b : core) {

        //R2_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, a, n, b);
        if (sj_ts != 0.0) {
          const auto Tnb = T.getv(n, b);
          const auto Sba = S.getv(b, a);
          const auto de_ts =
            (b.en() - n.en() + omega_t) * (a.en() - n.en() + omega_ts);
          Ak -= s_wa * sj_ts * Tnb * Sba * Wwavn / de_ts;
        }

        //R2_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, a, n, b);
        if (sj_st != 0.0) {
          const auto Snb = S.getv(n, b);
          const auto Tba = T.getv(b, a);
          const auto de_st =
            (b.en() - n.en() + omega_s) * (a.en() - n.en() + omega_ts);
          Ak -= s_wa * s_Kts * sj_st * Snb * Tba * Wwavn / de_st;
        }
      }
    }
  }
  return Ak / std::sqrt(2 * K + 1);
}

//==============================================================================

double dcp_L12(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
               double omega_s, int K, int kt, int ks,
               const Coulomb::meTable<double> &T,
               const Coulomb::meTable<double> &S,
               const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited) {

  // Both L1 and L2 share a Wk integral; faster to combine

  const auto s_Kts = Angular::neg1pow(K + kt + ks);
  const auto omega_ts = omega_t + omega_s;

  double Ak = 0.0;
#pragma omp parallel for reduction(+ : Ak) collapse(2)
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    for (std::size_t in = 0; in < excited.size(); ++in) {
      const auto &a = core[ia];
      const auto &n = excited[in];

      const auto s_wn = Angular::neg1pow_2(w.twoj() + n.twoj());
      const auto Wwnva = Coulomb::Wk_abcd(K, w, n, v, a);
      if (Wwnva == 0.0)
        continue;

      // L1:
      for (const auto &m : excited) {

        //L1_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, a, n, m);
        if (sj_ts != 0.0) {
          const auto Tmn = T.getv(m, n);
          const auto Sam = S.getv(a, m);
          const auto de_ts =
            (a.en() - n.en() - omega_ts) * (a.en() - m.en() - omega_s);
          Ak += s_wn * s_Kts * sj_ts * Tmn * Sam * Wwnva / de_ts;
        }

        //L1_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, a, n, m);
        if (sj_st != 0.0) {
          const auto Smn = S.getv(m, n);
          const auto Tam = T.getv(a, m);
          const auto de_st =
            (a.en() - n.en() - omega_ts) * (a.en() - m.en() - omega_t);
          Ak += s_wn * sj_st * Smn * Tam * Wwnva / de_st;
        }
      }

      // L2 (note -ve sign cf R1):
      for (const auto &b : core) {

        //L2_ts
        const auto sj_ts = Angular::SixJ(kt, ks, K, n, a, b);
        if (sj_ts != 0.0) {
          const auto Tab = T.getv(a, b);
          const auto Sbn = S.getv(b, n);
          const auto de_ts =
            (a.en() - n.en() - omega_ts) * (b.en() - n.en() - omega_s);
          Ak -= s_wn * sj_ts * Tab * Sbn * Wwnva / de_ts;
        }

        //L2_st
        const auto sj_st = ks == kt ? sj_ts : Angular::SixJ(ks, kt, K, n, a, b);
        if (sj_st != 0.0) {
          const auto Sab = S.getv(a, b);
          const auto Tbn = T.getv(b, n);
          const auto de_st =
            (a.en() - n.en() - omega_ts) * (b.en() - n.en() - omega_t);
          Ak -= s_wn * s_Kts * sj_st * Sab * Tbn * Wwnva / de_st;
        }
      }
    }
  }
  return Ak / std::sqrt(2 * K + 1);
}

//==============================================================================

double dcp(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
           double omega_s, int K, int kt, int ks,
           const Coulomb::meTable<double> &T, const Coulomb::meTable<double> &S,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited) {
  return dcp_C1(w, v, omega_t, omega_s, K, kt, ks, T, S, core, excited) +
         dcp_C2(w, v, omega_t, omega_s, K, kt, ks, T, S, core, excited) +
         dcp_R12(w, v, omega_t, omega_s, K, kt, ks, T, S, core, excited) +
         dcp_L12(w, v, omega_t, omega_s, K, kt, ks, T, S, core, excited);
}

} // namespace Amplitudes
