#include "Sigma2.hpp"

namespace MBPT {

//==============================================================================
double S_Sigma2(int k, const DiracSpinor &v, const DiracSpinor &w,
                const DiracSpinor &x, const DiracSpinor &y,
                const Coulomb::QkTable &qk,
                const std::vector<DiracSpinor> &core,
                const std::vector<DiracSpinor> &excited,
                const Angular::SixJTable *SJ = nullptr) {
  return S_Sigma2_a(k, v, w, x, y, qk, core, excited, SJ) +
         S_Sigma2_b(k, v, w, x, y, qk, core, excited, SJ) +
         S_Sigma2_c(k, v, w, x, y, qk, core, excited, SJ) +
         S_Sigma2_d(k, v, w, x, y, qk, core, excited, SJ);
}

//==============================================================================
double S_Sigma2_a(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable *SJ = nullptr) {
  /*
    S2^K_vwxy = [k]^-1 (-1)^k
                * (Qk_vnxa * Qk_awny + Pk_vnxa * Qk_awny + Qk_vnxa * Pk_awny)
                / (e_xa - e_vn)
  */

  const auto f = Angular::neg1pow(k) / (2.0 * k + 1.0);

  const auto de_xv = x.en() - v.en();

  double sum = 0.0;
  for (const auto &a : core) {
    for (const auto &n : excited) {
      const auto de = de_xv + a.en() - n.en();

      const auto qk_vnxa = qk.Q(k, v, n, x, a);
      const auto qk_awny = qk.Q(k, a, w, n, y);

      const auto pk_vnxa = qk.P(k, v, n, x, a);
      const auto pk_awny = qk.P(k, a, w, n, y);

      sum += (qk_vnxa * (qk_awny + pk_awny) + pk_vnxa * qk_awny) / de;
    }
  }

  return f * sum;
}

//==============================================================================
double S_Sigma2_b(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable *SJ = nullptr) {

  const auto f =
      Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj() + 2 * k) *
      (2.0 * k + 1.0);

  const auto de_yv = y.en() - v.en();

  if (!Coulomb::sixjTriads({}, {}, k, v, x, {}) ||
      !Coulomb::sixjTriads({}, {}, k, y, w, {}))
    return 0.0;

  double sum = 0.0;
  for (const auto &a : core) {
    for (const auto &n : excited) {

      const auto [u0, u1] = Coulomb::k_minmax_Q(v, n, a, y);
      const auto [l0, l1] = Coulomb::k_minmax_Q(a, w, x, n);
      if (l0 > l1)
        continue;

      if (!Coulomb::sixjTriads({}, {}, k, v, x, a) ||
          !Coulomb::sixjTriads({}, {}, k, y, w, n))
        continue;

      const auto de = de_yv + a.en() - n.en();

      for (int u = u0; u <= u1; u += 2) {
        const auto l0_sj = std::max(l0, std::abs(u - k));
        const auto l1_sj = std::min(l1, std::abs(u + k));
        for (int l = l0_sj; l <= l1_sj; l += 2) {

          // if (!Coulomb::triangle(l, u, k))
          //   continue;

          const auto sj1 = SJ->get(l, u, k, v, x, a);
          const auto sj2 = SJ->get(l, u, k, y, w, n);
          const auto s = Angular::neg1pow_2(2 * a.twoj() + 2 * l + 2 * u);

          const auto qk_vnay = qk.Q(u, v, n, a, y);
          const auto qk_awxn = qk.Q(l, a, w, x, n);

          sum += s * sj1 * sj2 * qk_vnay * qk_awxn / de;
        }
      }
    }
  }
  return f * sum;
}

//==============================================================================
double S_Sigma2_c(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable *SJ = nullptr) {

  const auto f =
      Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj() + 2 * k) *
      (2.0 * k + 1.0);

  const auto de_yv = x.en() - w.en();

  if (!Coulomb::sixjTriads({}, {}, k, v, x, {}) ||
      !Coulomb::sixjTriads({}, {}, k, y, w, {}))
    return 0.0;

  double sum = 0.0;
  for (const auto &a : core) {
    for (const auto &n : excited) {

      const auto [u0, u1] = Coulomb::k_minmax_Q(v, a, n, y);
      const auto [l0, l1] = Coulomb::k_minmax_Q(n, w, x, a);
      if (l0 > l1)
        continue;

      if (!Coulomb::sixjTriads({}, {}, k, v, x, n) ||
          !Coulomb::sixjTriads({}, {}, k, y, w, a))
        continue;

      const auto de = de_yv + a.en() - n.en();

      for (int u = u0; u <= u1; u += 2) {
        const auto l0_sj = std::max(l0, std::abs(u - k));
        const auto l1_sj = std::min(l1, std::abs(u + k));
        for (int l = l0_sj; l <= l1_sj; l += 2) {

          const auto sj1 = SJ->get(l, u, k, v, x, n);
          const auto sj2 = SJ->get(l, u, k, y, w, a);
          const auto s = Angular::neg1pow_2(2 * a.twoj() + 2 * l + 2 * u);

          const auto qk_vany = qk.Q(u, v, a, n, y);
          const auto qk_nwxa = qk.Q(l, n, w, x, a);

          sum += s * sj1 * sj2 * qk_vany * qk_nwxa / de;
        }
      }
    }
  }
  return f * sum;
}

//==============================================================================
double S_Sigma2_d(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable *SJ = nullptr) {

  const auto f =
      Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj() + 2 * k) *
      (2.0 * k + 1.0);

  const auto de_vw = -v.en() - w.en();

  if (!Coulomb::sixjTriads({}, {}, k, v, x, {}) ||
      !Coulomb::sixjTriads({}, {}, k, w, y, {})) // check!
    return 0.0;

  double sum = 0.0;
  for (const auto &a : core) {
    for (const auto &b : core) {

      const auto [u0, u1] = Coulomb::k_minmax_Q(v, w, a, b);
      const auto [l0, l1] = Coulomb::k_minmax_Q(a, b, x, y);
      if (l0 > l1)
        continue;

      if (!Coulomb::sixjTriads({}, {}, k, v, x, a) ||
          !Coulomb::sixjTriads({}, {}, k, w, y, b))
        continue;

      const auto de = de_vw + a.en() + b.en();

      for (int u = u0; u <= u1; u += 2) {
        const auto l0_sj = std::max(l0, std::abs(u - k));
        const auto l1_sj = std::min(l1, std::abs(u + k));
        for (int l = l0_sj; l <= l1_sj; l += 2) {

          const auto sj1 = SJ->get(l, u, k, v, x, a);
          const auto sj2 = SJ->get(l, u, k, w, y, b);
          const auto s = Angular::neg1pow_2(2 * a.twoj() + 2 * l + 2 * u);

          const auto qk_vwab = qk.Q(u, v, w, a, b);
          const auto qk_abxy = qk.Q(l, a, b, x, y);

          sum += s * sj1 * sj2 * qk_vwab * qk_abxy / de;
        }
      }
    }
  }
  return f * sum;
}

} // namespace MBPT