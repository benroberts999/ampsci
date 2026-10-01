#include "Sigma2.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cassert>

namespace MBPT {

//==============================================================================
// Returns string representation of Denominators enum
std::string parse_Denominators(Denominators d) {
  switch (d) {
  case Denominators::RS:
    return "RS";
  case Denominators::Fermi:
    return "Fermi";
  case Denominators::Fermi0:
    return "Fermi0";
  case Denominators::DFK:
    return "DFK";
  case Denominators::BW:
    return "BW";
  }
  assert(false);
  return "Error";
}

// Parses string to Denominators enum (case-insensitive); returns DFK if unrecognised
Denominators parse_Denominators(std::string_view s) {
  if (qip::ci_compare(s, "RS"))
    return Denominators::RS;
  if (qip::ci_compare(s, "Fermi"))
    return Denominators::Fermi;
  if (qip::ci_compare(s, "Fermi0"))
    return Denominators::Fermi0;
  if (qip::ci_compare(s, "BW"))
    return Denominators::BW;
  return Denominators::DFK;
}

//==============================================================================
// External-leg part of a Sigma_2 energy denominator; see Sigma2.hpp for the
// full description of the modes and parameters. Full denominator is
// (e_a - e_n) + leg_de. For BW, the whole denominator is E0 - E_intermediate,
// where E_intermediate = es + ei + e_n - e_a is the total energy of the state
// between the two Coulomb vertices; hence the target slot becomes (E0 - es).
double leg_de(Denominators denominators, double et_bar, double et,
              double ei_bar, double ei, double es, double E0) {
  switch (denominators) {
  case Denominators::RS:
    return et - ei;
  case Denominators::Fermi:
    return et_bar - ei_bar;
  case Denominators::Fermi0:
    return 0.0;
  case Denominators::DFK:
    return et_bar - ei;
  case Denominators::BW:
    return (E0 - es) - ei;
  }
  assert(false);
  return 0.0;
}

//==============================================================================
double e_bar(int kappa_v, const std::vector<DiracSpinor> &excited) {
  // Assumes excited is sorted by energy, so first match of each kappa is the
  // lowest (true in practice).
  const auto v_bar = std::find_if(
    excited.cbegin(), excited.cend(),
    [kappa_v](const DiracSpinor &n) { return n.kappa() == kappa_v; });
  if (v_bar == excited.cend()) {
    return 0.0;
  }
  return v_bar->en();
}

//==============================================================================
bool Sk_vwxy_SR(int k, const DiracSpinor &v, const DiracSpinor &w,
                const DiracSpinor &x, const DiracSpinor &y) {
  return Coulomb::sixjTriads({}, {}, k, v, x, {}) &&
         Coulomb::sixjTriads({}, {}, k, w, y, {});
}

//==============================================================================
std::pair<int, int> k_minmax_S(const DiracSpinor &v, const DiracSpinor &w,
                               const DiracSpinor &x, const DiracSpinor &y) {

  return k_minmax_S(v.twoj(), w.twoj(), x.twoj(), y.twoj());
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
               const Angular::SixJTable &SixJ, Denominators denominators,
               const std::vector<double> &fk, double E0) {
  using namespace Sigma2;

  if (!Sk_vwxy_SR(k, v, w, x, y))
    return 0.0;

  return S_Sigma2_ab(k, v, w, x, y, qk, core, excited, SixJ, denominators, fk,
                     E0) +
         S_Sigma2_c1(k, v, w, x, y, qk, core, excited, SixJ, denominators, fk,
                     E0) +
         S_Sigma2_c2(k, v, w, x, y, qk, core, excited, SixJ, denominators, fk,
                     E0) +
         S_Sigma2_d(k, v, w, x, y, qk, core, excited, SixJ, denominators, fk,
                    E0);
}

//==============================================================================
//==============================================================================
//==============================================================================
double Sigma2::S_Sigma2_ab(int k, const DiracSpinor &v, const DiracSpinor &w,
                           const DiracSpinor &x, const DiracSpinor &y,
                           const Coulomb::QkTable &qk,
                           const std::vector<DiracSpinor> &core,
                           const std::vector<DiracSpinor> &excited,
                           const Angular::SixJTable &SixJ,
                           Denominators denominators,
                           const std::vector<double> &fk, double E0) {

  // overall selectrion rule tested outside

  // screening factors
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  const auto f = Angular::neg1pow(k) / (2.0 * k + 1.0);

  const auto v0 = e_bar(v.kappa(), excited);
  const auto w0 = e_bar(w.kappa(), excited);
  const auto x0 = e_bar(x.kappa(), excited);
  const auto y0 = e_bar(y.kappa(), excited);

  // External-leg parts of the energy denominators (full denominator is
  // (e_a - e_n) + de_ext); _hc is that of the bra-ket partner:
  // diagram a: (e_x - e_v) [target x, intermediate v]; partner: (e_w - e_y)
  // diagram b: (e_y - e_w); partner: (e_v - e_x)
  // The intermediate state of diagram a is {v, y, n, a}, and of diagram b is
  // {x, w, n, a}: the es argument is the valence orbital that is not the
  // intermediate leg (only used for BW)
  const auto de_ext_a =
    leg_de(denominators, x0, x.en(), v0, v.en(), y.en(), E0);
  const auto de_ext_a_hc =
    leg_de(denominators, w0, w.en(), y0, y.en(), v.en(), E0);
  const auto de_ext_b =
    leg_de(denominators, y0, y.en(), w0, w.en(), x.en(), E0);
  const auto de_ext_b_hc =
    leg_de(denominators, v0, v.en(), x0, x.en(), w.en(), E0);

  double sum = 0.0;
  for (const auto &a : core) {
    for (const auto &n : excited) {
      const auto de_an = a.en() - n.en();

      // Hermitianse: average each diagram with its dagged partner
      // [diagram a of <xy|S|vw> has the numerator of diagram b of <vw|S|xy>,
      // with the external part of the denominator negated, and vice versa].
      // Gives symmetric CI matrix; correct to this order.
      const auto inv_de_a =
        0.5 * (1.0 / (de_an + de_ext_a) + 1.0 / (de_an + de_ext_a_hc));
      const auto inv_de_b =
        0.5 * (1.0 / (de_an + de_ext_b) + 1.0 / (de_an + de_ext_b_hc));

      // A diagrams:
      const auto qk_vnxa = qk.Q(k, v, n, x, a);
      const auto pk_vnxa = qk.P2(k, v, n, x, a, SixJ, fk);

      const auto qk_awny = qk.Q(k, a, w, n, y);
      const auto pk_awny = qk.P2(k, a, w, n, y, SixJ, fk);
      const auto wk_awny = qk_awny + pk_awny;

      // Screening: the direct leg is of multipolarity k, so each term picks up
      // a single Fk(k); the exchange legs are screened inside P2 (at their own
      // internal multipolarity). The bubble (screening) diagrams a1 and b1 are
      // thus screened once, all others twice.

      // diagrams a1, a2, a3:
      sum += Fk(k) * (qk_vnxa * wk_awny + pk_vnxa * qk_awny) * inv_de_a;

      // B diagrams: a <-> n
      const auto qk_vaxn = qk_vnxa;
      const auto pk_vaxn = v == x ? pk_vnxa : qk.P2(k, v, a, x, n, SixJ, fk);
      const auto qk_nway = qk_awny;
      const auto pk_nway = w == y ? pk_awny : qk.P2(k, n, w, a, y, SixJ, fk);
      const auto wk_nway = qk_nway + pk_nway;

      // diagrams b1, b2, b3:
      sum += Fk(k) * (qk_vaxn * wk_nway + pk_vaxn * qk_nway) * inv_de_b;
    }
  }

  return f * sum;
}

//==============================================================================
double Sigma2::S_Sigma2_c1(int k, const DiracSpinor &v, const DiracSpinor &w,
                           const DiracSpinor &x, const DiracSpinor &y,
                           const Coulomb::QkTable &qk,
                           const std::vector<DiracSpinor> &core,
                           const std::vector<DiracSpinor> &excited,
                           const Angular::SixJTable &SixJ,
                           Denominators denominators,
                           const std::vector<double> &fk, double E0) {

  // overall selectrion rule tested outside

  // screening factors
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  const auto f =
    Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj() + 2 * k) *
    (2.0 * k + 1.0);

  const auto v0 = e_bar(v.kappa(), excited);
  const auto w0 = e_bar(w.kappa(), excited);
  const auto x0 = e_bar(x.kappa(), excited);
  const auto y0 = e_bar(y.kappa(), excited);

  // External-leg parts of the energy denominators (full denominator is
  // (e_a - e_n) + de_ext):
  // c1's own denominator has (e_y - e_v) [target y, intermediate v];
  // its bra-ket partner (c1 of the reversed element, same numerator)
  // has (e_w - e_x). The intermediate state is {x, v, n, a}
  const auto de_ext = leg_de(denominators, y0, y.en(), v0, v.en(), x.en(), E0);
  const auto de_ext_hc =
    leg_de(denominators, w0, w.en(), x0, x.en(), v.en(), E0);

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

      const auto de_an = a.en() - n.en();
      // Hermitise: average with bra-ket partner (see S_Sigma2_ab)
      const auto inv_de =
        0.5 * (1.0 / (de_an + de_ext) + 1.0 / (de_an + de_ext_hc));

      for (int u = u0; u <= u1; u += 2) {
        const auto l0_SixJ = l0; // allow += 2
        const auto l1_SixJ = std::min(l1, std::abs(u + k));
        for (int l = l0_SixJ; l <= l1_SixJ; l += 2) {

          const auto SixJ1 = SixJ.get(l, u, k, v, x, a);
          const auto SixJ2 = SixJ.get(l, u, k, y, w, n);
          const auto s = Angular::neg1pow_2(2 * a.twoj() + 2 * l + 2 * u);

          const auto qk_vnay = Fk(u) * qk.Q(u, v, n, a, y);
          const auto qk_awxn = Fk(l) * qk.Q(l, a, w, x, n);

          sum += s * SixJ1 * SixJ2 * qk_vnay * qk_awxn * inv_de;
        }
      }
    }
  }
  return f * sum;
}

//==============================================================================
double Sigma2::S_Sigma2_c2(int k, const DiracSpinor &v, const DiracSpinor &w,
                           const DiracSpinor &x, const DiracSpinor &y,
                           const Coulomb::QkTable &qk,
                           const std::vector<DiracSpinor> &core,
                           const std::vector<DiracSpinor> &excited,
                           const Angular::SixJTable &SixJ,
                           Denominators denominators,
                           const std::vector<double> &fk, double E0) {

  // overall selectrion rule tested outside

  // screening factors
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  const auto v0 = e_bar(v.kappa(), excited);
  const auto w0 = e_bar(w.kappa(), excited);
  const auto x0 = e_bar(x.kappa(), excited);
  const auto y0 = e_bar(y.kappa(), excited);

  const auto f =
    Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj() + 2 * k) *
    (2.0 * k + 1.0);

  // External-leg parts of the energy denominators (full denominator is
  // (e_a - e_n) + de_ext):
  // c2's denominator has (e_x - e_w) [target x, intermediate w];
  // its HC partner has (e_v - e_y). The intermediate state is {w, y, n, a}
  const auto de_ext = leg_de(denominators, x0, x.en(), w0, w.en(), y.en(), E0);
  const auto de_ext_hc =
    leg_de(denominators, v0, v.en(), y0, y.en(), w.en(), E0);

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

      const auto de_an = a.en() - n.en();
      // Hermitianise: average with HC partner
      const auto inv_de =
        0.5 * (1.0 / (de_an + de_ext) + 1.0 / (de_an + de_ext_hc));

      for (int u = u0; u <= u1; u += 2) {
        const auto l0_SixJ = l0; // allow += 2
        const auto l1_SixJ = std::min(l1, std::abs(u + k));
        for (int l = l0_SixJ; l <= l1_SixJ; l += 2) {

          const auto SixJ1 = SixJ.get(l, u, k, v, x, n);
          const auto SixJ2 = SixJ.get(l, u, k, y, w, a);
          const auto s = Angular::neg1pow_2(2 * a.twoj() + 2 * l + 2 * u);

          const auto qk_vany = Fk(u) * qk.Q(u, v, a, n, y);
          const auto qk_nwxa = Fk(l) * qk.Q(l, n, w, x, a);

          sum += s * SixJ1 * SixJ2 * qk_vany * qk_nwxa * inv_de;
        }
      }
    }
  }
  return f * sum;
}

//==============================================================================
double Sigma2::S_Sigma2_d(int k, const DiracSpinor &v, const DiracSpinor &w,
                          const DiracSpinor &x, const DiracSpinor &y,
                          const Coulomb::QkTable &qk,
                          const std::vector<DiracSpinor> &core,
                          const std::vector<DiracSpinor> &excited,
                          const Angular::SixJTable &SixJ,
                          Denominators denominators,
                          const std::vector<double> &fk, double E0) {

  // screening factors
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  const auto f = Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj()) *
                 (2.0 * k + 1.0);

  const auto v0 = e_bar(v.kappa(), excited);
  const auto w0 = e_bar(w.kappa(), excited);
  const auto x0 = e_bar(x.kappa(), excited);
  const auto y0 = e_bar(y.kappa(), excited);
  const auto e0 = DiracSpinor::min_En(excited);

  // The intermediate state of diagram d holds all four valence orbitals
  // (plus the two holes), so BW gives (e_v + e_w + e_x + e_y) - E0
  const auto e_vwxy = v.en() + w.en() + x.en() + y.en();

  // Here, "Fermi0" cancellation doesn't happen: valence energies remain.
  // Own denominator has -(e_v + e_w); bra-ket partner has -(e_x + e_y)
  // Both external legs are target-state legs, so DFK coincides with Fermi.
  const auto pair_en = [denominators, e0, e_vwxy, E0](
                         double ea_bar, double ea, double eb_bar, double eb) {
    switch (denominators) {
    case Denominators::RS:
      return ea + eb;
    case Denominators::Fermi:
    case Denominators::DFK:
      return ea_bar + eb_bar;
    case Denominators::Fermi0:
      return 2.0 * e0;
    case Denominators::BW:
      return e_vwxy - E0;
    }
    assert(false);
    return 0.0;
  };
  const auto e_vw = pair_en(v0, v.en(), w0, w.en());
  const auto e_xy = pair_en(x0, x.en(), y0, y.en());

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

      const auto e_ab = a.en() + b.en();
      // Hermitianise: average with bra-ket partner (see S_Sigma2_ab)
      const auto inv_de = 0.5 * (1.0 / (e_ab - e_vw) + 1.0 / (e_ab - e_xy));
      const auto s = Angular::neg1pow_2(a.twoj() - b.twoj());

      for (int u = u0; u <= u1; u += 2) {
        const auto l0_SixJ = l0;
        const auto l1_SixJ = std::min(l1, std::abs(u + k));
        for (int l = l0_SixJ; l <= l1_SixJ; l += 2) {

          if (!Coulomb::triangle(l, u, k))
            continue;

          const auto SixJ1 = SixJ.get(l, u, k, v, x, a);
          const auto SixJ2 = SixJ.get(l, u, k, w, y, b);

          const auto qu_vwab = Fk(u) * qk.Q(u, v, w, a, b);
          const auto ql_abxy = Fk(l) * qk.Q(l, a, b, x, y);

          sum += s * SixJ1 * SixJ2 * qu_vwab * ql_abxy * inv_de;
        }
      }
    }
  }
  return f * sum;
}

//==============================================================================
Coulomb::LkTable calculate_Sk(const std::string &filename,
                              const std::vector<DiracSpinor> &external,
                              const std::vector<DiracSpinor> &core,
                              const std::vector<DiracSpinor> &excited,
                              const Coulomb::QkTable &qk, int max_k,
                              bool exclude_wrong_parity_box,
                              Denominators denominators, bool no_new_integrals,
                              const std::vector<double> &fk, double E0) {

  Coulomb::LkTable Sk;

  const auto max_twoj =
    std::max({DiracSpinor::max_tj(excited), DiracSpinor::max_tj(core),
              DiracSpinor::max_tj(external)});
  Angular::SixJTable sjt(max_twoj);

  const auto Sk_function = [&](int k, const DiracSpinor &v,
                               const DiracSpinor &w, const DiracSpinor &x,
                               const DiracSpinor &y) {
    return MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited, sjt, denominators,
                         fk, E0);
  };
  const auto Sk_selection_rule = [&](int k, const DiracSpinor &v,
                                     const DiracSpinor &w, const DiracSpinor &x,
                                     const DiracSpinor &y) {
    return exclude_wrong_parity_box ? Coulomb::Qk_abcd_SR(k, v, w, x, y) :
                                      MBPT::Sk_vwxy_SR(k, v, w, x, y);
  };

  // Try to read from disk (may already have calculated Qk)
  Sk.read(filename);

  const auto existing = Sk.count();
  {
    if (!no_new_integrals)
      Sk.fill(external, Sk_function, Sk_selection_rule, max_k);

    const auto total = Sk.count();
    assert(total >= existing);
    const auto new_integrals = total - existing;
    std::cout << "Calculated " << new_integrals << " new MBPT integrals\n";
    if (new_integrals > 0) {
      Sk.write(filename);
    }
  }
  std::cout << "\n" << std::flush;
  return Sk;
}

//==============================================================================
std::vector<double> average_hk(const Coulomb::LkTable &Sk,
                               const Coulomb::QkTable &qk,
                               const std::vector<DiracSpinor> &external,
                               int max_k) {
  // Only integrals with |Q^k| above this contribute to the average:
  const double q_cut = 1.0e-3;

  const auto tj_max = DiracSpinor::max_tj(external);
  const auto k_max = max_k < 0 ? tj_max : std::min(max_k, tj_max);

  std::vector<double> hk(std::size_t(k_max + 1), 0.0);
  std::vector<double> sum(std::size_t(k_max + 1), 0.0);
  std::vector<long> num(std::size_t(k_max + 1), 0);

  // S^k has a 4-fold symmetry: S^k_vwxy = S^k_wvyx = S^k_xyvw = S^k_yxwv.
  // The Lk normal order accounts for the first pair; the bra-ket pair is added
  // here. Count each distinct integral once - otherwise each is weighted by its
  // own symmetry multiplicity, which is not the same for every integral.
  const auto normal_order = [&Sk](const DiracSpinor &v, const DiracSpinor &w,
                                  const DiracSpinor &x, const DiracSpinor &y) {
    return std::min(Sk.NormalOrder(v, w, x, y), Sk.NormalOrder(x, y, v, w));
  };

  // Parallel over k: each k accumulates independently
#pragma omp parallel for
  for (int k = 0; k <= k_max; ++k) {
    for (const auto &v : external) {
      for (const auto &w : external) {
        for (const auto &x : external) {
          for (const auto &y : external) {
            const auto index = normal_order(v, w, x, y);
            if (index != Sk.CurrentOrder(v, w, x, y))
              continue;
            const auto [k0, k1] = Coulomb::k_minmax_Q(v, w, x, y);
            if (k < k0 || k > k1 || (k - k0) % 2 != 0)
              continue;
            const auto s = Sk.Q(k, index);
            if (s == 0.0)
              continue;
            const auto q = qk.Q(k, v, w, x, y);
            if (std::abs(q) < q_cut)
              continue;
            sum[std::size_t(k)] += s / q;
            ++num[std::size_t(k)];
          }
        }
      }
    }
  }
  for (std::size_t k = 0; k < hk.size(); ++k) {
    if (num[k] > 0) {
      hk[k] = sum[k] / double(num[k]);
    }
  }
  return hk;
}

} // namespace MBPT
