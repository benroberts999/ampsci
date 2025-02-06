#include "Sigma2.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"

namespace MBPT {

//==============================================================================
std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>>
split_basis(const std::vector<DiracSpinor> &basis, double E_Fermi,
            int min_n_core, int max_n_excited) {

  std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>> core_excited;
  auto &[core, excited] = core_excited;
  for (const auto &Fn : basis) {
    if (Fn.en() <= E_Fermi && Fn.n() >= min_n_core) {
      core.push_back(Fn);
    } else if (Fn.en() > E_Fermi && Fn.n() <= max_n_excited) {
      excited.push_back(Fn);
    }
  }
  return core_excited;
}

//==============================================================================
double e_bar(int kappa_v, const std::vector<DiracSpinor> &excited) {
  const auto v_bar = std::find_if(
      excited.cbegin(), excited.cend(),
      [kappa_v](const DiracSpinor &n) { return n.kappa() == kappa_v; });
  assert(v_bar != excited.cend());
  return v_bar->en();
}

//==============================================================================
bool Sk_vwxy_SR(int k, const DiracSpinor &v, const DiracSpinor &w,
                const DiracSpinor &x, const DiracSpinor &y) {
  return Coulomb::sixjTriads({}, {}, k, v, x, {}) &&
         Coulomb::sixjTriads({}, {}, k, w, y, {});
}

//==============================================================================
int number_below_Fermi(const DiracSpinor &i, const DiracSpinor &j,
                       const DiracSpinor &k, const DiracSpinor &l,
                       double eFermi) {
  int num_core = 0;
  if (i.en() < eFermi)
    ++num_core;
  if (j.en() < eFermi)
    ++num_core;
  if (k.en() < eFermi)
    ++num_core;
  if (l.en() < eFermi)
    ++num_core;
  return num_core;
}

//==============================================================================
std::pair<int, int> k_minmax_S(const DiracSpinor &v, const DiracSpinor &w,
                               const DiracSpinor &x, const DiracSpinor &y) {

  // From the 6j part only:
  // |b-d| <= k <=|b+d|
  // |a-c| <= k <=|a+c|
  const auto [lk1, uk1] = Coulomb::k_minmax_tj(v.twoj(), x.twoj());
  const auto [lk2, uk2] = Coulomb::k_minmax_tj(w.twoj(), y.twoj());

  return {std::max({lk1, lk2}), std::min({uk1, uk2})};
}

std::pair<int, int> k_minmax_S(int twojv, int twojw, int twojx, int twojy) {

  // From the 6j part only:
  // |b-d| <= k <=|b+d|
  // |a-c| <= k <=|a+c|
  const auto [lk1, uk1] = Coulomb::k_minmax_tj(twojv, twojx);
  const auto [lk2, uk2] = Coulomb::k_minmax_tj(twojw, twojy);

  return {std::max({lk1, lk2}), std::min({uk1, uk2})};
}

//==============================================================================
double Sk_vwxy(int k, const DiracSpinor &v, const DiracSpinor &w,
               const DiracSpinor &x, const DiracSpinor &y,
               const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited,
               const Angular::SixJTable &SixJ, Denominators denominators) {
  using namespace InternalSigma;

  if (!Sk_vwxy_SR(k, v, w, x, y))
    return 0.0;

  return S_Sigma2_ab(k, v, w, x, y, qk, core, excited, SixJ, denominators) +
         S_Sigma2_c1(k, v, w, x, y, qk, core, excited, SixJ, denominators) +
         S_Sigma2_c2(k, v, w, x, y, qk, core, excited, SixJ, denominators) +
         S_Sigma2_d(k, v, w, x, y, qk, core, excited, SixJ, denominators);
}

//==============================================================================
//==============================================================================
//==============================================================================
double InternalSigma::S_Sigma2_ab(int k, const DiracSpinor &v,
                                  const DiracSpinor &w, const DiracSpinor &x,
                                  const DiracSpinor &y,
                                  const Coulomb::QkTable &qk,
                                  const std::vector<DiracSpinor> &core,
                                  const std::vector<DiracSpinor> &excited,
                                  const Angular::SixJTable &SixJ,
                                  Denominators denominators) {

  // overall selectrion rule tested outside

  const auto f = Angular::neg1pow(k) / (2.0 * k + 1.0);

  // 1. Use actual RS denoms
  // 2. Use symmetrised RS denoms
  // 3. Use "lowest kappa" denoms
  // 4. Symmetrised "lowest kappa" denoms
  // 5. "Full BW approx": only na

  // const auto v0 = e_bar(v.kappa(), excited);
  // const auto w0 = e_bar(w.kappa(), excited);
  // const auto x0 = e_bar(x.kappa(), excited);
  // const auto y0 = e_bar(y.kappa(), excited);

  // const auto de_xv = x.en() - v.en();
  const auto de_xv = denominators == Denominators::BW ?
                         0.0 :
                         //  0.5 * (x0 - v0 + y0 - w0) :
                         0.5 * (x.en() - v.en() + y.en() - w.en());

  double sum = 0.0;
  for (const auto &a : core) {
    for (const auto &n : excited) {
      const auto de = de_xv + a.en() - n.en();

      // A diagrams:
      const auto qk_vnxa = qk.Q(k, v, n, x, a);
      const auto pk_vnxa = qk.P(k, v, n, x, a, &SixJ);

      const auto qk_awny = qk.Q(k, a, w, n, y);
      const auto pk_awny = qk.P(k, a, w, n, y, &SixJ);
      const auto wk_awny = qk_awny + pk_awny;

      // diagrams a1, a2, a3:
      sum += (qk_vnxa * wk_awny + pk_vnxa * qk_awny) / de;

      // B diagrams: a <-> n
      const auto qk_vaxn = qk_vnxa;
      const auto pk_vaxn = v == x ? pk_vnxa : qk.P(k, v, a, x, n, &SixJ);
      const auto qk_nway = qk_awny;
      const auto pk_nway = w == y ? pk_awny : qk.P(k, n, w, a, y, &SixJ);
      const auto wk_nway = qk_nway + pk_nway;

      // diagrams b1, b2, b3:
      sum += (qk_vaxn * wk_nway + pk_vaxn * qk_nway) / de;
    }
  }

  return f * sum;
}

//==============================================================================
double InternalSigma::S_Sigma2_c1(int k, const DiracSpinor &v,
                                  const DiracSpinor &w, const DiracSpinor &x,
                                  const DiracSpinor &y,
                                  const Coulomb::QkTable &qk,
                                  const std::vector<DiracSpinor> &core,
                                  const std::vector<DiracSpinor> &excited,
                                  const Angular::SixJTable &SixJ,
                                  Denominators denominators) {

  // overall selectrion rule tested outside

  const auto f =
      Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj() + 2 * k) *
      (2.0 * k + 1.0);

  // const auto v0 = e_bar(v.kappa(), excited);
  // const auto w0 = e_bar(w.kappa(), excited);
  // const auto x0 = e_bar(x.kappa(), excited);
  // const auto y0 = e_bar(y.kappa(), excited);

  // const auto de_yv = y.en() - v.en();
  const auto de_yv = denominators == Denominators::BW ?
                         //  0.5 * (y0 - v0 + x0 - w0) :
                         0.0 :
                         0.5 * (y.en() - v.en() + x.en() - w.en());

  double sum = 0.0;
  for (const auto &a : core) {
    if (!Coulomb::sixjTriads({}, {}, k, v, x, a))
      continue;
    for (const auto &n : excited) {
      const auto [u0, u1] = Coulomb::k_minmax_Q(v, n, a, y);
      const auto [l0, l1] = Coulomb::k_minmax_Q(a, w, x, n);
      if (l0 > l1)
        continue;

      if (!Coulomb::sixjTriads({}, {}, k, y, w, n))
        continue;

      const auto de = de_yv + a.en() - n.en();

      for (int u = u0; u <= u1; u += 2) {
        const auto l0_SixJ = l0; // allow += 2
        const auto l1_SixJ = std::min(l1, std::abs(u + k));
        for (int l = l0_SixJ; l <= l1_SixJ; l += 2) {

          const auto SixJ1 = SixJ.get(l, u, k, v, x, a);
          const auto SixJ2 = SixJ.get(l, u, k, y, w, n);
          const auto s = Angular::neg1pow_2(2 * a.twoj() + 2 * l + 2 * u);

          const auto qk_vnay = qk.Q(u, v, n, a, y);
          const auto qk_awxn = qk.Q(l, a, w, x, n);

          sum += s * SixJ1 * SixJ2 * qk_vnay * qk_awxn / de;
        }
      }
    }
  }
  return f * sum;
}

//==============================================================================
double InternalSigma::S_Sigma2_c2(int k, const DiracSpinor &v,
                                  const DiracSpinor &w, const DiracSpinor &x,
                                  const DiracSpinor &y,
                                  const Coulomb::QkTable &qk,
                                  const std::vector<DiracSpinor> &core,
                                  const std::vector<DiracSpinor> &excited,
                                  const Angular::SixJTable &SixJ,
                                  Denominators denominators) {

  // overall selectrion rule tested outside

  // const auto v0 = e_bar(v.kappa(), excited);
  // const auto w0 = e_bar(w.kappa(), excited);
  // const auto x0 = e_bar(x.kappa(), excited);
  // const auto y0 = e_bar(y.kappa(), excited);

  const auto f =
      Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj() + 2 * k) *
      (2.0 * k + 1.0);

  const auto de_yv = denominators == Denominators::BW ?
                         0.0 :
                         //  0.5 * (x0 - w0 + y0 - v0) :
                         0.5 * (x.en() - w.en() + y.en() - v.en());

  double sum = 0.0;
  for (const auto &a : core) {
    if (!Coulomb::sixjTriads({}, {}, k, y, w, a))
      continue;
    for (const auto &n : excited) {

      const auto [u0, u1] = Coulomb::k_minmax_Q(v, a, n, y);
      const auto [l0, l1] = Coulomb::k_minmax_Q(n, w, x, a);
      if (l0 > l1)
        continue;

      if (!Coulomb::sixjTriads({}, {}, k, v, x, n))
        continue;

      const auto de = de_yv + a.en() - n.en();

      for (int u = u0; u <= u1; u += 2) {
        const auto l0_SixJ = l0; // allow += 2
        const auto l1_SixJ = std::min(l1, std::abs(u + k));
        for (int l = l0_SixJ; l <= l1_SixJ; l += 2) {

          const auto SixJ1 = SixJ.get(l, u, k, v, x, n);
          const auto SixJ2 = SixJ.get(l, u, k, y, w, a);
          const auto s = Angular::neg1pow_2(2 * a.twoj() + 2 * l + 2 * u);

          const auto qk_vany = qk.Q(u, v, a, n, y);
          const auto qk_nwxa = qk.Q(l, n, w, x, a);

          sum += s * SixJ1 * SixJ2 * qk_vany * qk_nwxa / de;
        }
      }
    }
  }
  return f * sum;
}

//==============================================================================
double InternalSigma::S_Sigma2_d(int k, const DiracSpinor &v,
                                 const DiracSpinor &w, const DiracSpinor &x,
                                 const DiracSpinor &y,
                                 const Coulomb::QkTable &qk,
                                 const std::vector<DiracSpinor> &core,
                                 const std::vector<DiracSpinor> &excited,
                                 const Angular::SixJTable &SixJ,
                                 Denominators denominators) {

  const auto f = Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj()) *
                 (2.0 * k + 1.0);

  const auto vbar = e_bar(v.kappa(), excited);
  const auto wbar = e_bar(w.kappa(), excited);
  const auto xbar = e_bar(x.kappa(), excited);
  const auto ybar = e_bar(y.kappa(), excited);

  // symmetrised..
  const auto de_vw = denominators == Denominators::BW ?
                         -0.5 * (vbar + wbar + xbar + ybar) :
                         -0.5 * (v.en() + w.en() + x.en() + y.en());

  double sum = 0.0;
  for (const auto &a : core) {
    if (!Coulomb::sixjTriads({}, {}, k, v, x, a))
      continue;

    for (const auto &b : core) {

      const auto [u0, u1] = Coulomb::k_minmax_Q(v, w, a, b);
      const auto [l0, l1] = Coulomb::k_minmax_Q(a, b, x, y);
      if (l0 > l1)
        continue;

      if (!Coulomb::sixjTriads({}, {}, k, w, y, b))
        continue;

      const auto de = de_vw + a.en() + b.en();
      const auto s = Angular::neg1pow_2(a.twoj() - b.twoj());

      for (int u = u0; u <= u1; u += 2) {
        const auto l0_SixJ = l0;
        const auto l1_SixJ = std::min(l1, std::abs(u + k));
        for (int l = l0_SixJ; l <= l1_SixJ; l += 2) {

          if (!Coulomb::triangle(l, u, k))
            continue;

          const auto SixJ1 = SixJ.get(l, u, k, v, x, a);
          const auto SixJ2 = SixJ.get(l, u, k, w, y, b);

          const auto qu_vwab = qk.Q(u, v, w, a, b);
          const auto ql_abxy = qk.Q(l, a, b, x, y);

          sum += s * SixJ1 * SixJ2 * qu_vwab * ql_abxy / de;
        }
      }
    }
  }
  return f * sum;
}

} // namespace MBPT