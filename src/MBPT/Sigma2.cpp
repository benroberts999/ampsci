#include "Sigma2.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "MBPT/Feynman.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <optional>

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
  if (qip::ci_compare(s, "DFK"))
    return Denominators::DFK;
  std::cout << "Unknown denominator option: " << s << "? Defaulting to DFK\n;";
  return Denominators::DFK;
}

//==============================================================================
// External-leg part of a Sigma_2 energy denominator; see Sigma2.hpp for the
// full description of the modes and parameters. Full denominator is
// (e_a - e_n) + leg_de. For BW, the whole denominator is E0 - E_intermediate,
// where E_intermediate = es + ei + e_n - e_a is the total energy of the state
// between the two Coulomb vertices; hence the target slot becomes (E0 - es).
CIOrbitals::CIOrbitals(const std::vector<DiracSpinor> &holes,
                       const std::vector<DiracSpinor> &ci_basis) {
  for (const auto &orbitals : {&holes, &ci_basis}) {
    for (const auto &orbital : *orbitals) {
      const auto index = std::size_t(orbital.nk_index());
      if (index >= m_member.size()) {
        m_member.resize(index + 1, false);
      }
      m_member[index] = true;
    }
  }
}

//==============================================================================
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
  // Each vertex conserves parity, so S^k vanishes unless the pair states vw
  // and xy have the same parity; k itself may have either parity
  return v.parity() * w.parity() == x.parity() * y.parity() &&
         Coulomb::sixjTriads({}, {}, k, v, x, {}) &&
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
double e_bar(const DiracSpinor &leg, const std::vector<DiracSpinor> &core,
             const std::vector<DiracSpinor> &excited) {
  const auto hole_leg =
    std::find(core.cbegin(), core.cend(), leg) != core.cend();
  return hole_leg ? leg.en() : e_bar(leg.kappa(), excited);
}

//==============================================================================
double Sk_vwxy(int k, const DiracSpinor &v, const DiracSpinor &w,
               const DiracSpinor &x, const DiracSpinor &y,
               const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited,
               const Angular::SixJTable &SixJ, Denominators denominators,
               const std::vector<double> &fk, double E0,
               const std::vector<RMatrix> *screening_dQk,
               const CIOrbitals *ci_orbitals) {
  using namespace Sigma2;

  if (!Sk_vwxy_SR(k, v, w, x, y))
    return 0.0;

  // Hole legs: external legs that are core orbitals (the holes of hole CSFs)
  const auto is_hole_leg = [&core](const DiracSpinor &leg) {
    return std::find(core.cbegin(), core.cend(), leg) != core.cend();
  };
  const bool hole_leg_vw = is_hole_leg(v) || is_hole_leg(w);
  const bool hole_leg_xy = is_hole_leg(x) || is_hole_leg(y);
  // The bra and ket themselves are intermediate states of terms with a hole
  // leg (zero denominators); they lie within the CI orbitals
  assert((!(hole_leg_vw || hole_leg_xy) || ci_orbitals != nullptr) &&
         "Sigma_2 with a hole leg requires the CI orbitals");

  // Static loop diagrams from the Feynman method, for the k the correction
  // covers; otherwise all diagrams are second-order Goldstone. Not with the CI
  // orbitals left out: the loop's intermediate states cannot be restricted
  const auto ik = std::size_t(k);
  const bool feynman_loop = screening_dQk != nullptr &&
                            ik < screening_dQk->size() &&
                            ci_orbitals == nullptr;

  const auto Sk_Goldstone = [&](const DiracSpinor &v_, const DiracSpinor &w_,
                                const DiracSpinor &x_, const DiracSpinor &y_) {
    return S_Sigma2_ab(k, v_, w_, x_, y_, qk, core, excited, SixJ, denominators,
                       fk, E0, feynman_loop, ci_orbitals) +
           S_Sigma2_c1(k, v_, w_, x_, y_, qk, core, excited, SixJ, denominators,
                       fk, E0, ci_orbitals) +
           S_Sigma2_c2(k, v_, w_, x_, y_, qk, core, excited, SixJ, denominators,
                       fk, E0, ci_orbitals) +
           S_Sigma2_d(k, v_, w_, x_, y_, qk, core, excited, SixJ, denominators,
                      fk, E0, ci_orbitals) +
           // The particle-particle ladder, with a hole leg in each pair
           (hole_leg_vw && hole_leg_xy ?
              S_Sigma2_pp(k, v_, w_, x_, y_, qk, core, excited, SixJ,
                          denominators, E0, ci_orbitals) :
              0.0);
  };

  // With the Fermi-level denominators, a hole leg is one more orbital of its
  // kappa: at the energy of the lowest excited orbital of that kappa, in
  // every slot (as for its Sigma_1)
  const bool holes_at_kappa_energy =
    (hole_leg_vw || hole_leg_xy) &&
    (denominators == Denominators::DFK || denominators == Denominators::Fermi);
  if (holes_at_kappa_energy) {
    const auto at_kappa_energy = [&](const DiracSpinor &leg) {
      auto out = leg;
      const auto lowest =
        std::find_if(excited.cbegin(), excited.cend(),
                     [kappa = leg.kappa()](const DiracSpinor &n) {
                       return n.kappa() == kappa;
                     });
      if (is_hole_leg(leg) && lowest != excited.cend()) {
        out.en() = lowest->en();
      }
      return out;
    };
    return Sk_Goldstone(at_kappa_energy(v), at_kappa_energy(w),
                        at_kappa_energy(x), at_kappa_energy(y));
  }

  return feynman_loop ?
           Sk_Goldstone(v, w, x, y) +
             S_Sigma2_screen(k, v, w, x, y, screening_dQk->at(ik)) :
           Sk_Goldstone(v, w, x, y);
}

//==============================================================================
//==============================================================================
//==============================================================================
double Sigma2::S_Sigma2_ab(
  int k, const DiracSpinor &v, const DiracSpinor &w, const DiracSpinor &x,
  const DiracSpinor &y, const Coulomb::QkTable &qk,
  const std::vector<DiracSpinor> &core, const std::vector<DiracSpinor> &excited,
  const Angular::SixJTable &SixJ, Denominators denominators,
  const std::vector<double> &fk, double E0, bool exclude_static_loop,
  const CIOrbitals *ci_orbitals) {

  // overall selectrion rule tested outside

  // Every term has a direct Coulomb line of multipolarity k between a core and
  // an excited orbital, Q^k_vnxa or Q^k_awny, whose other end is (v, x) or
  // (w, y). If k is not allowed at either of those (k of the 'wrong' parity,
  // which only the box diagrams have), every term vanishes
  if (!Angular::Ck_kk_SR(k, v.kappa(), x.kappa()) &&
      !Angular::Ck_kk_SR(k, w.kappa(), y.kappa()))
    return 0.0;

  // screening factors
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  const auto f = Angular::neg1pow(k) / (2.0 * k + 1.0);

  const auto v0 = e_bar(v, core, excited);
  const auto w0 = e_bar(w, core, excited);
  const auto x0 = e_bar(x, core, excited);
  const auto y0 = e_bar(y, core, excited);

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
      // ...and that line vanishes unless k is allowed between a and n
      if (!Angular::Ck_kk_SR(k, a.kappa(), n.kappa()))
        continue;

      // Both internal lines within the CI orbitals: left out
      if (ci_orbitals && ci_orbitals->contains(a) && ci_orbitals->contains(n))
        continue;

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
      // internal multipolarity). The screening diagrams a1 and b1 are
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

      // The loop screening diagrams a1 and b1 are the Q*Q parts of the above.
      // With
      // exclude_static_loop, their static part (denominator e_a - e_n,
      // i.e. 1/de_an, for each) is removed, with the same screening factor:
      // only the energy dependence of their denominators remains, screened by
      // fk as the other diagrams are. The static loop, screened to all
      // orders, is added by the caller
      if (exclude_static_loop) {
        const auto qq_loop = qk_vnxa * qk_awny;
        sum -= Fk(k) * qq_loop * 2.0 / de_an;
      }
    }
  }

  return f * sum;
}

//==============================================================================
double Sigma2::S_Sigma2_c1(
  int k, const DiracSpinor &v, const DiracSpinor &w, const DiracSpinor &x,
  const DiracSpinor &y, const Coulomb::QkTable &qk,
  const std::vector<DiracSpinor> &core, const std::vector<DiracSpinor> &excited,
  const Angular::SixJTable &SixJ, Denominators denominators,
  const std::vector<double> &fk, double E0, const CIOrbitals *ci_orbitals) {

  // overall selectrion rule tested outside

  // screening factors
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  const auto f =
    Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj() + 2 * k) *
    (2.0 * k + 1.0);

  const auto v0 = e_bar(v, core, excited);
  const auto w0 = e_bar(w, core, excited);
  const auto x0 = e_bar(x, core, excited);
  const auto y0 = e_bar(y, core, excited);

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
      // Both internal lines within the CI orbitals: left out
      if (ci_orbitals && ci_orbitals->contains(a) && ci_orbitals->contains(n))
        continue;
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
          const auto s =
            Angular::neg1pow_2(a.twoj() + n.twoj() + 2 * l + 2 * u);

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
double Sigma2::S_Sigma2_c2(
  int k, const DiracSpinor &v, const DiracSpinor &w, const DiracSpinor &x,
  const DiracSpinor &y, const Coulomb::QkTable &qk,
  const std::vector<DiracSpinor> &core, const std::vector<DiracSpinor> &excited,
  const Angular::SixJTable &SixJ, Denominators denominators,
  const std::vector<double> &fk, double E0, const CIOrbitals *ci_orbitals) {

  // overall selectrion rule tested outside

  // screening factors
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  const auto v0 = e_bar(v, core, excited);
  const auto w0 = e_bar(w, core, excited);
  const auto x0 = e_bar(x, core, excited);
  const auto y0 = e_bar(y, core, excited);

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
      // Both internal lines within the CI orbitals: left out
      if (ci_orbitals && ci_orbitals->contains(a) && ci_orbitals->contains(n))
        continue;

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
          const auto s =
            Angular::neg1pow_2(a.twoj() + n.twoj() + 2 * l + 2 * u);

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
double Sigma2::S_Sigma2_d(
  int k, const DiracSpinor &v, const DiracSpinor &w, const DiracSpinor &x,
  const DiracSpinor &y, const Coulomb::QkTable &qk,
  const std::vector<DiracSpinor> &core, const std::vector<DiracSpinor> &excited,
  const Angular::SixJTable &SixJ, Denominators denominators,
  const std::vector<double> &fk, double E0, const CIOrbitals *ci_orbitals) {

  // screening factors
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  const auto f = Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj()) *
                 (2.0 * k + 1.0);

  const auto v0 = e_bar(v, core, excited);
  const auto w0 = e_bar(w, core, excited);
  const auto x0 = e_bar(x, core, excited);
  const auto y0 = e_bar(y, core, excited);
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
      // Both internal holes within the CI orbitals: left out
      if (ci_orbitals && ci_orbitals->contains(a) && ci_orbitals->contains(b))
        continue;

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
double Sigma2::S_Sigma2_pp(int k, const DiracSpinor &v, const DiracSpinor &w,
                           const DiracSpinor &x, const DiracSpinor &y,
                           const Coulomb::QkTable &qk,
                           const std::vector<DiracSpinor> &core,
                           const std::vector<DiracSpinor> &excited,
                           const Angular::SixJTable &SixJ,
                           Denominators denominators, double E0,
                           const CIOrbitals *ci_orbitals) {

  const auto f = Angular::neg1pow_2(v.twoj() + w.twoj() + x.twoj() + y.twoj()) *
                 (2.0 * k + 1.0);

  const auto v0 = e_bar(v, core, excited);
  const auto w0 = e_bar(w, core, excited);
  const auto x0 = e_bar(x, core, excited);
  const auto y0 = e_bar(y, core, excited);
  const auto e0 = DiracSpinor::min_En(excited);

  // Both pairs are consumed at the vertices (target-state legs); the
  // intermediate state holds only the two excited orbitals
  const auto pair_en = [denominators, e0, E0](double ea_bar, double ea,
                                              double eb_bar, double eb) {
    switch (denominators) {
    case Denominators::RS:
      return ea + eb;
    case Denominators::Fermi:
    case Denominators::DFK:
      return ea_bar + eb_bar;
    case Denominators::Fermi0:
      return 2.0 * e0;
    case Denominators::BW:
      return E0;
    }
    assert(false);
    return 0.0;
  };
  const auto e_vw = pair_en(v0, v.en(), w0, w.en());
  const auto e_xy = pair_en(x0, x.en(), y0, y.en());

  double sum = 0.0;
  for (const auto &m : excited) {
    if (!Coulomb::sixjTriads({}, {}, k, v, x, m))
      continue;

    for (const auto &n : excited) {
      // Both internal lines within the CI orbitals: left out
      if (ci_orbitals && ci_orbitals->contains(m) && ci_orbitals->contains(n))
        continue;

      const auto [u0, u1] = Coulomb::k_minmax_Q(v, w, m, n);
      const auto [l0, l1] = Coulomb::k_minmax_Q(m, n, x, y);
      if (l0 > l1)
        continue;

      if (!Coulomb::sixjTriads({}, {}, k, w, y, n))
        continue;

      const auto e_mn = m.en() + n.en();
      // Hermitianise: average with bra-ket partner (see S_Sigma2_ab)
      const auto inv_de = 0.5 * (1.0 / (e_xy - e_mn) + 1.0 / (e_vw - e_mn));
      const auto s = Angular::neg1pow_2(m.twoj() - n.twoj());

      for (int u = u0; u <= u1; u += 2) {
        const auto l1_SixJ = std::min(l1, std::abs(u + k));
        for (int l = l0; l <= l1_SixJ; l += 2) {

          if (!Coulomb::triangle(l, u, k))
            continue;

          const auto SixJ1 = SixJ.get(l, u, k, v, x, m);
          const auto SixJ2 = SixJ.get(l, u, k, w, y, n);

          const auto qu_vwmn = qk.Q(u, v, w, m, n);
          const auto ql_mnxy = qk.Q(l, m, n, x, y);

          sum += s * SixJ1 * SixJ2 * qu_vwmn * ql_mnxy * inv_de;
        }
      }
    }
  }
  return f * sum;
}

//==============================================================================
double Sigma2::S_Sigma2_screen(int k, const DiracSpinor &v,
                               const DiracSpinor &w, const DiracSpinor &x,
                               const DiracSpinor &y, const RMatrix &delta_Qk) {

  // Same angular structure as Q^k_vwxy (the loop diagrams have a direct
  // Coulomb line of multipolarity k at each vertex), with the radial Coulomb
  // kernel replaced by the screening kernel
  if (!Angular::Ck_kk_SR(k, v.kappa(), x.kappa()) ||
      !Angular::Ck_kk_SR(k, w.kappa(), y.kappa()))
    return 0.0;

  const auto angular = Angular::neg1pow_2(2 * k + v.twoj() - w.twoj()) *
                       Angular::Ck_kk(k, v.kappa(), x.kappa()) *
                       Angular::Ck_kk(k, w.kappa(), y.kappa());

  return angular * two_body_ME(delta_Qk, v, w, x, y);
}

//==============================================================================
Coulomb::WkTable calculate_Sk(
  const std::string &filename, const std::vector<DiracSpinor> &external,
  const std::vector<DiracSpinor> &pair_basis,
  const std::vector<DiracSpinor> &core, const std::vector<DiracSpinor> &excited,
  const Coulomb::QkTable &qk, int max_k, bool exclude_wrong_parity_box,
  Denominators denominators, bool no_new_integrals,
  const std::vector<double> &fk, double E0,
  const std::vector<RMatrix> *screening_dQk, const CIOrbitals *ci_orbitals,
  const InCISpace *hole_CSFs) {

  Coulomb::WkTable Sk;

  const auto max_twoj =
    std::max({DiracSpinor::max_tj(excited), DiracSpinor::max_tj(core),
              DiracSpinor::max_tj(external)});
  Angular::SixJTable sjt(max_twoj);

  const auto Sk_function = [&](int k, const DiracSpinor &v,
                               const DiracSpinor &w, const DiracSpinor &x,
                               const DiracSpinor &y) {
    return MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited, sjt, denominators,
                         fk, E0, screening_dQk, ci_orbitals);
  };

  // Which states are in pair_basis, and which are core (hole legs), indexed
  // by nk_index: checked for every candidate integral, so must be fast
  const auto membership = [](const std::vector<DiracSpinor> &set) {
    std::vector<bool> in_set;
    for (const auto &a : set) {
      const auto index = std::size_t(a.nk_index());
      if (index >= in_set.size()) {
        in_set.resize(index + 1, false);
      }
      in_set[index] = true;
    }
    return in_set;
  };
  const auto in_pair_basis = membership(pair_basis);
  const auto in_core = membership(core);
  const auto is_in = [](const std::vector<bool> &in_set, const DiracSpinor &a) {
    const auto ia = std::size_t(a.nk_index());
    return ia < in_set.size() && in_set[ia];
  };
  const auto pair_in_basis = [&](const DiracSpinor &a, const DiracSpinor &b) {
    return is_in(in_pair_basis, a) && is_in(in_pair_basis, b);
  };

  // Either pair, (v, w) or (x, y), in pair_basis; or diagonal (xy in the
  // same shells as vw), which enter the energy of every configuration.
  // Diagonal is by non-relativistic shell (n, l), both j: otherwise the j
  // components of one shell would be treated differently, making Sigma_2 spin
  // dependent (it would then split LS multiplets). Unchanged under the
  // symmetries of S^k (vwxy = wvyx = xyvw = yxwv)
  const auto same_shell = [](const DiracSpinor &a, const DiracSpinor &b) {
    return a.n() == b.n() && a.l() == b.l();
  };
  // Hole legs (core orbitals among the external legs): the integrals the
  // block's hole CSFs need (hole_CSFs). One hole leg: the other pair are the
  // two particles of its hole CSF. Two hole legs (one in each pair; two in
  // one pair never occurs, as the pairs of the CSF matrix elements are
  // particle-particle or particle-hole): each hole with the particle of the
  // other pair in a hole CSF (the spectator is not known). The bra and ket
  // states are then in the CI space, so the terms with a zero denominator
  // are among those left out. Nothing without the predicate. Unchanged
  // under the symmetries of S^k (which swap within pairs, or the pairs).
  // Returns nothing if there is no hole leg
  const auto hole_leg_rule = [&](const DiracSpinor &v, const DiracSpinor &w,
                                 const DiracSpinor &x,
                                 const DiracSpinor &y) -> std::optional<bool> {
    const auto hole_v = is_in(in_core, v);
    const auto hole_w = is_in(in_core, w);
    const auto hole_x = is_in(in_core, x);
    const auto hole_y = is_in(in_core, y);
    const auto holes_vw = int(hole_v) + int(hole_w);
    const auto holes_xy = int(hole_x) + int(hole_y);
    if (holes_vw + holes_xy == 0)
      return std::nullopt;
    if (holes_vw > 1 || holes_xy > 1 || hole_CSFs == nullptr)
      return false;
    if (holes_vw == 1 && holes_xy == 1) {
      const auto &hole_in_vw = hole_v ? v : w;
      const auto &particle_in_vw = hole_v ? w : v;
      const auto &hole_in_xy = hole_x ? x : y;
      const auto &particle_in_xy = hole_x ? y : x;
      return (*hole_CSFs)(hole_in_xy, particle_in_vw, nullptr) &&
             (*hole_CSFs)(hole_in_vw, particle_in_xy, nullptr);
    }
    if (holes_vw == 1)
      return (*hole_CSFs)(hole_v ? v : w, x, &y);
    return (*hole_CSFs)(hole_x ? x : y, v, &w);
  };

  const auto Sk_selection_rule = [&](int k, const DiracSpinor &v,
                                     const DiracSpinor &w, const DiracSpinor &x,
                                     const DiracSpinor &y) {
    const auto hole_legs = hole_leg_rule(v, w, x, y);
    if (hole_legs) {
      if (!*hole_legs)
        return false;
    } else {
      const auto diagonal = (same_shell(x, v) && same_shell(y, w)) ||
                            (same_shell(x, w) && same_shell(y, v));
      if (!diagonal && !pair_in_basis(v, w) && !pair_in_basis(x, y))
        return false;
    }
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
    std::cout << "Calculated " << new_integrals << " new MBPT integrals\n"
              << std::flush;
    if (new_integrals > 0) {
      Sk.write(filename);
    }
  }
  std::cout << "\n" << std::flush;
  return Sk;
}

//==============================================================================
std::vector<double> average_hk(const Coulomb::WkTable &Sk,
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

  // S^k has a 4-fold symmetry: S^k_vwxy = S^k_wvyx = S^k_xyvw = S^k_yxwv,
  // which the table's normal order accounts for. Count each distinct integral
  // once - otherwise each is weighted by its own symmetry multiplicity, which
  // is not the same for every integral.

  // Parallel over k: each k accumulates independently
#pragma omp parallel for
  for (int k = 0; k <= k_max; ++k) {
    for (const auto &v : external) {
      for (const auto &w : external) {
        for (const auto &x : external) {
          for (const auto &y : external) {
            const auto index = Sk.NormalOrder(v, w, x, y);
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
