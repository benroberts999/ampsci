#include "CI_Integrals.hpp"
#include "CSF.hpp"
#include "Coulomb/include.hpp"
#include "HF/Breit.hpp"
#include "LinAlg/Matrix.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/Sigma2.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/format.hpp"
#include "fmt/ostream.hpp"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace CI {

//==============================================================================
// Calculates the anti-symmetrised Coulomb integral for 2-particle states:
// C1*C2*(g_abcd-g_abdc), where Cs are C.G. coefficients
double CSF2_Coulomb(const Coulomb::QkTable &qk, DiracSpinor::Index v,
                    DiracSpinor::Index w, DiracSpinor::Index x,
                    DiracSpinor::Index y, int twoJ) {

  // If c==d, or a==b : can make short-cut due to symmetry
  // More efficient to use two Q's than W:

  const auto tjv = Angular::nkindex_to_twoj(v);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tjx = Angular::nkindex_to_twoj(x);
  const auto tjy = Angular::nkindex_to_twoj(y);

  const auto kv = Angular::nkindex_to_kappa(v);
  const auto kw = Angular::nkindex_to_kappa(w);
  const auto kx = Angular::nkindex_to_kappa(x);
  const auto ky = Angular::nkindex_to_kappa(y);

  double out = 0.0;

  // Direct part:
  const auto [k0, k1] = Coulomb::k_minmax_Q(kv, kw, kx, ky);
  for (int k = k0; k <= k1; k += 2) {
    const auto sjs = Angular::sixj_2(tjv, tjw, twoJ, tjy, tjx, 2 * k);
    if (sjs == 0.0)
      continue;
    const auto qk_abcd = qk.Q(k, v, w, x, y);
    const auto s = Angular::neg1pow_2(tjv + tjx + 2 * k + twoJ);
    out += s * sjs * qk_abcd;
  }
  // const auto f = Angular::neg1pow_2(tjv + tjx + 2 * k0 + twoJ) / (twoJ + 1.0);
  // out += f * qk.P(twoJ / 2, v, y, w, x);

  // Take advantage of symmetries: faster (+ numerically stable)
  // c == d => J is even (identical states), eta2=1/sqrt(2)
  // eta_ab = 1/sqrt(2) if a==b
  // Therefore: e.g., if c==d
  // => eta_ab * eta_cd * (out + (-1)^J*out) = eta_ab * sqrt(2) * out
  if (v == w && x == y) {
    return out;
  } else if (x == y || v == w) {
    // by {ab},{cd} symmetry: same works for case a==b
    return std::sqrt(2.0) * out;
  }

  // Exchange part:
  const auto [l0, l1] = Coulomb::k_minmax_Q(kv, kw, ky, kx);
  for (int k = l0; k <= l1; k += 2) {
    const auto sjs = Angular::sixj_2(tjv, tjw, twoJ, tjx, tjy, 2 * k);
    if (sjs == 0.0)
      continue;
    const auto qk_abdc = qk.Q(k, v, w, y, x);
    const auto s = Angular::neg1pow_2(tjv + tjx + 2 * k);
    out += s * sjs * qk_abdc;
  }

  return out;
}

//==============================================================================
double CSF2_Sigma2(const Coulomb::WkTable &Sk, DiracSpinor::Index v,
                   DiracSpinor::Index w, DiracSpinor::Index x,
                   DiracSpinor::Index y, int twoJ, const Coulomb::QkTable *qk,
                   const std::vector<double> &hk, const Coulomb::WkTable *dSk,
                   double dE0) {

  // If c==d, or a==b : can make short-cut due to symmetry
  // More efficient to use two Q's than W:

  double out = 0.0;

  const auto tjv = Angular::nkindex_to_twoj(v);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tjx = Angular::nkindex_to_twoj(x);
  const auto tjy = Angular::nkindex_to_twoj(y);

  // S^k for a diagram that was not calculated (i.e., outside the cis2 basis):
  // extrapolate as S^k = hk*Q^k (rather than storing the extrapolated table).
  // These get no E0 shift: there is no derivative to shift them with.
  const auto extrap = [qk, &hk](int k, DiracSpinor::Index a,
                                DiracSpinor::Index b, DiracSpinor::Index c,
                                DiracSpinor::Index d) {
    if (qk == nullptr || std::size_t(k) >= hk.size())
      return 0.0;
    return hk[std::size_t(k)] * qk->Q(k, a, b, c, d);
  };

  // S^k, shifted to the target-level E0 if the derivative table is given
  const auto get_Sk =
    [&Sk, dSk, dE0, &extrap](int k, DiracSpinor::Index a, DiracSpinor::Index b,
                             DiracSpinor::Index c, DiracSpinor::Index d) {
      if (!Sk.contains(k, a, b, c, d)) {
        return extrap(k, a, b, c, d);
      }
      const auto S0 = Sk.Q(k, a, b, c, d);
      if (dSk == nullptr) {
        return S0;
      }
      return corrected_Sk(S0, dSk->Q(k, a, b, c, d), dE0);
    };

  // Direct part:
  const auto [k0, k1] = MBPT::k_minmax_S(tjv, tjw, tjx, tjy);
  for (int k = k0; k <= k1; ++k) {
    const auto sjs = Angular::sixj_2(tjv, tjw, twoJ, tjy, tjx, 2 * k);
    if (sjs == 0.0)
      continue;
    const auto Sk_abcd = get_Sk(k, v, w, x, y);
    const auto s = Angular::neg1pow_2(tjv + tjx + 2 * k + twoJ);
    out += s * sjs * Sk_abcd;
  }

  // Take advantage of symmetries: faster (+ numerically stable)
  // c == d => J is even (identical states), eta2=1/sqrt(2)
  // eta_ab = 1/sqrt(2) if a==b
  // Therefore: e.g., if c==d
  // => eta_ab * eta_cd * (out + (-1)^J*out) = eta_ab * sqrt(2) * out
  if (v == w && x == y) {
    return out;
  } else if (x == y || v == w) {
    // by {ab},{cd} symmetry: same works for case a==b
    return std::sqrt(2.0) * out;
  }

  // Exchange part:
  const auto [l0, l1] = MBPT::k_minmax_S(tjv, tjw, tjy, tjx);
  for (int k = l0; k <= l1; ++k) {
    const auto sjs = Angular::sixj_2(tjv, tjw, twoJ, tjx, tjy, 2 * k);
    if (sjs == 0.0)
      continue;
    const auto Sk_abdc = get_Sk(k, v, w, y, x);
    const auto s = Angular::neg1pow_2(tjv + tjx + 2 * k);
    out += s * sjs * Sk_abdc;
  }

  return out;
}

//==============================================================================
// Calculates the anti-symmetrised Breit integral for 2-particle states:
// C1*C2*(b_abcd-b_abdc), where Cs are C.G. coefficients
double CSF2_Breit(const Coulomb::WkTable &Bk, DiracSpinor::Index v,
                  DiracSpinor::Index w, DiracSpinor::Index x,
                  DiracSpinor::Index y, int twoJ) {

  // If c==d, or a==b : can make short-cut due to symmetry
  // More efficient to use two Q's than W:

  const auto tjv = Angular::nkindex_to_twoj(v);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tjx = Angular::nkindex_to_twoj(x);
  const auto tjy = Angular::nkindex_to_twoj(y);

  double out = 0.0;

  // Direct part:
  const auto [k0, k1] = HF::Breit::k_minmax_tj(tjv, tjw, tjx, tjy);
  for (int k = k0; k <= k1; ++k) {
    const auto sjs = Angular::sixj_2(tjv, tjw, twoJ, tjy, tjx, 2 * k);
    if (sjs == 0.0)
      continue;
    const auto bk_abcd = Bk.Q(k, v, w, x, y);
    const auto s = Angular::neg1pow_2(tjv + tjx + 2 * k + twoJ);
    out += s * sjs * bk_abcd;
  }

  // Take advantage of symmetries: faster (+ numerically stable)
  if (v == w && x == y) {
    return out;
  } else if (x == y || v == w) {
    // by {ab},{cd} symmetry: same works for case a==b
    return std::sqrt(2.0) * out;
  }

  // Exchange part:
  const auto [l0, l1] = HF::Breit::k_minmax_tj(tjv, tjw, tjy, tjx);
  for (int k = l0; k <= l1; ++k) {
    const auto sjs = Angular::sixj_2(tjv, tjw, twoJ, tjx, tjy, 2 * k);
    if (sjs == 0.0)
      continue;
    const auto bk_abdc = Bk.Q(k, v, w, y, x);
    const auto s = Angular::neg1pow_2(tjv + tjx + 2 * k);
    out += s * sjs * bk_abdc;
  }

  return out;
}

//==============================================================================
double Sigma2_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                 const Coulomb::WkTable &Sk, const Coulomb::QkTable *qk,
                 const std::vector<double> &hk, const Coulomb::WkTable *dSk,
                 double dE0) {
  const auto [v, w] = A.states;
  const auto [x, y] = B.states;
  return CSF2_Sigma2(Sk, v, w, x, y, twoJ, qk, hk, dSk, dE0);
}

//==============================================================================
double Breit_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                const Coulomb::WkTable &Bk) {
  const auto [v, w] = A.states;
  const auto [x, y] = B.states;
  return CSF2_Breit(Bk, v, w, x, y, twoJ);
}

//==============================================================================
double Sigma2_AB(const CI::CSF3 &A, const CI::CSF3 &B, int,
                 const Coulomb::WkTable &Sk, const Coulomb::QkTable *qk,
                 const std::vector<double> &hk, const Coulomb::WkTable *dSk,
                 double dE0,
                 const std::map<DiracSpinor::Index, double> *orbital_energies) {
  assert((dSk == nullptr || orbital_energies != nullptr) &&
         "Brillouin-Wigner Sigma_2 requires the orbital energies");

  if (CSF3::num_different(A, B) > 2)
    return 0.0;

  // Two-body: acts on electrons 1 and 2; the third electron must coincide
  double out = 0.0;
  for (const auto &p : A.parents()) {
    for (const auto &q : B.parents()) {
      if (p.spectator != q.spectator || p.twoJ2 != q.twoJ2)
        continue;
      // The pair energy excludes the third electron
      const auto dE0_pair =
        dSk == nullptr ? dE0 : dE0 - orbital_energies->at(p.spectator);
      out += p.coef * q.coef *
             CSF2_Sigma2(Sk, p.pair[0], p.pair[1], q.pair[0], q.pair[1],
                         p.twoJ2, qk, hk, dSk, dE0_pair);
    }
  }
  return 3.0 * out;
}

//==============================================================================
double Breit_AB(const CI::CSF3 &A, const CI::CSF3 &B, int,
                const Coulomb::WkTable &Bk) {
  if (CSF3::num_different(A, B) > 2)
    return 0.0;

  // Two-body: acts on electrons 1 and 2; the third electron must coincide
  double out = 0.0;
  for (const auto &p : A.parents()) {
    for (const auto &q : B.parents()) {
      if (p.spectator != q.spectator || p.twoJ2 != q.twoJ2)
        continue;
      out +=
        p.coef * q.coef *
        CSF2_Breit(Bk, p.pair[0], p.pair[1], q.pair[0], q.pair[1], p.twoJ2);
    }
  }
  return 3.0 * out;
}

//==============================================================================
double Sigma2_AB(const CI::CSF2p1h &A, const CI::CSF2p1h &B, int twoJ,
                 const Coulomb::WkTable &Sk, const Coulomb::QkTable *qk,
                 const std::vector<double> &hk, const Coulomb::WkTable *dSk,
                 double dE0,
                 const std::map<DiracSpinor::Index, double> *orbital_energies) {
  assert((dSk == nullptr || orbital_energies != nullptr) &&
         "Brillouin-Wigner Sigma_2 requires the orbital energies");
  // The pair energy excludes the hole, whose energy enters with the opposite
  // sign
  const auto dE0_pair =
    dSk == nullptr ? dE0 : dE0 + orbital_energies->at(B.hole);
  return two_body_2p1h(sigma2_pair_element(Sk, qk, hk, dSk, dE0_pair), A, B,
                       twoJ);
}

//==============================================================================
double Sigma2_AB(const CI::CSF1 &A, const CI::CSF2p1h &B, int twoJ,
                 const Coulomb::WkTable &Sk, const Coulomb::QkTable *qk,
                 const std::vector<double> &hk, const Coulomb::WkTable *dSk,
                 double dE0,
                 const std::map<DiracSpinor::Index, double> *orbital_energies) {
  assert((dSk == nullptr || orbital_energies != nullptr) &&
         "Brillouin-Wigner Sigma_2 requires the orbital energies");
  const auto dE0_pair =
    dSk == nullptr ? dE0 : dE0 + orbital_energies->at(B.hole);
  return two_body_1p_2p1h(sigma2_pair_element(Sk, qk, hk, dSk, dE0_pair), A, B,
                          twoJ);
}

//==============================================================================
double Breit_AB(const CI::CSF2p1h &A, const CI::CSF2p1h &B, int twoJ,
                const Coulomb::WkTable &Bk) {
  return two_body_2p1h(breit_pair_element(Bk), A, B, twoJ);
}

//==============================================================================
double Breit_AB(const CI::CSF1 &A, const CI::CSF2p1h &B, int twoJ,
                const Coulomb::WkTable &Bk) {
  return two_body_1p_2p1h(breit_pair_element(Bk), A, B, twoJ);
}

//==============================================================================
double Sigma2_AB(const CI::CSF &A, const CI::CSF &B, int twoJ,
                 const Coulomb::WkTable &Sk, const Coulomb::QkTable *qk,
                 const std::vector<double> &hk, const Coulomb::WkTable *dSk,
                 double dE0,
                 const std::map<DiracSpinor::Index, double> *orbital_energies) {
  return std::visit(
    [&](const auto &a, const auto &b) -> double {
      using A_type = std::decay_t<decltype(a)>;
      using B_type = std::decay_t<decltype(b)>;
      if constexpr (std::is_same_v<A_type, CSF1> &&
                    std::is_same_v<B_type, CSF1>) {
        // No valence pairs
        return 0.0;
      } else if constexpr (std::is_same_v<A_type, CSF2> &&
                           std::is_same_v<B_type, CSF2>) {
        return Sigma2_AB(a, b, twoJ, Sk, qk, hk, dSk, dE0);
      } else if constexpr ((std::is_same_v<A_type, CSF3> &&
                            std::is_same_v<B_type, CSF3>) ||
                           (std::is_same_v<A_type, CSF2p1h> &&
                            std::is_same_v<B_type, CSF2p1h>) ||
                           (std::is_same_v<A_type, CSF1> &&
                            std::is_same_v<B_type, CSF2p1h>)) {
        return Sigma2_AB(a, b, twoJ, Sk, qk, hk, dSk, dE0, orbital_energies);
      } else if constexpr (std::is_same_v<A_type, CSF2p1h> &&
                           std::is_same_v<B_type, CSF1>) {
        // Transposed element (the CI matrix is symmetric)
        return Sigma2_AB(b, a, twoJ, Sk, qk, hk, dSk, dE0, orbital_energies);
      } else {
        assert(false && "Sigma_2 between CSFs of different kinds is not "
                        "implemented");
        return 0.0;
      }
    },
    A.variant(), B.variant());
}

//==============================================================================
double Breit_AB(const CI::CSF &A, const CI::CSF &B, int twoJ,
                const Coulomb::WkTable &Bk) {
  return std::visit(
    [&](const auto &a, const auto &b) -> double {
      using A_type = std::decay_t<decltype(a)>;
      using B_type = std::decay_t<decltype(b)>;
      if constexpr (std::is_same_v<A_type, CSF1> &&
                    std::is_same_v<B_type, CSF1>) {
        // No valence pairs
        return 0.0;
      } else if constexpr (std::is_same_v<A_type, B_type> ||
                           (std::is_same_v<A_type, CSF1> &&
                            std::is_same_v<B_type, CSF2p1h>)) {
        return Breit_AB(a, b, twoJ, Bk);
      } else if constexpr (std::is_same_v<A_type, CSF2p1h> &&
                           std::is_same_v<B_type, CSF1>) {
        // Transposed element (the CI matrix is symmetric)
        return Breit_AB(b, a, twoJ, Bk);
      } else {
        assert(false && "Breit between CSFs of different kinds is not "
                        "implemented");
        return 0.0;
      }
    },
    A.variant(), B.variant());
}

//==============================================================================
double corrected_Sigma(double Sigma, double dSigma, double dE) {
  if (Sigma == 0.0) {
    return 0.0;
  }
  const auto denom = Sigma - dE * dSigma;
  if (denom == 0.0) {
    return Sigma;
  }
  const auto Sigma_c = Sigma * Sigma / denom;
  // The resummed (Pade) form has a pole at dE = Sigma/dSigma. If the
  // correction enhances |Sigma|, dE*dSigma is approaching Sigma (the pole),
  // where the formula diverges: distrust, return uncorrected Sigma
  return std::abs(Sigma_c) > std::abs(Sigma) ? Sigma : Sigma_c;
}

//==============================================================================
double corrected_Sk(double Sk, double dSk, double dE0) {
  if (Sk == 0.0) {
    return 0.0;
  }
  const auto denom = Sk - dE0 * dSk;
  // The resummed (Pade) form has a single pole, at dE0 = Sk/dSk. For a
  // single-denominator term it is the physical degeneracy (the intermediate
  // state crossing E0); for the actual multi-term sum, once dE0 crosses the
  // nearest pole the one-pole form no longer approximates anything. A sign
  // change of the denominator means we are past it: return unshifted value
  if (denom * Sk <= 0.0) {
    return Sk;
  }
  return Sk * Sk / denom;
}

//==============================================================================
double
Sigma1Correction::delta_h1(DiracSpinor::Index a, DiracSpinor::Index b,
                           std::optional<DiracSpinor::Index> spectator1,
                           std::optional<DiracSpinor::Index> spectator2) const {
  const auto Sig = S1.getv(a, b);
  // E0 = 0.0 means it was never set: no correction (see iterate_E0)
  if (Sig == 0.0 || E0 == 0.0) {
    return 0.0;
  }
  const auto es = e_sigma.find(Angular::nkindex_to_kappa(a));
  if (es == e_sigma.end()) {
    return 0.0;
  }
  auto dE = E0 - es->second;
  for (const auto &spectator : {spectator1, spectator2}) {
    if (!spectator)
      continue;
    const auto e_spec = en.find(*spectator);
    if (e_spec == en.end())
      return 0.0;
    dE -= e_spec->second;
  }
  return corrected_Sigma(Sig, dS1.getv(a, b), dE) - Sig;
}

//==============================================================================
bool Sigma1Correction::read_write(const std::string &fname, IO::FRW::RoW rw) {
  using IO::FRW::rw_binary;
  const bool readQ = (rw == IO::FRW::read);

  if (readQ && !IO::FRW::file_exists(fname)) {
    return false;
  }

  std::fstream f;
  IO::FRW::open_binary(f, fname, rw);
  if (!f) {
    return false;
  }

  // Each container: size, then (key, value) pairs.
  // nb: E0 is not stored: it belongs to a single (J, parity) and is set at
  // solve time (see iterate_E0)
  const auto rw_map = [&f, rw, readQ](auto *map) {
    using Map = std::decay_t<decltype(*map)>;
    auto size = map->size();
    rw_binary(f, rw, size);
    if (readQ) {
      map->clear();
      for (std::size_t i = 0; i < size; ++i) {
        typename Map::key_type key{};
        typename Map::mapped_type value{};
        rw_binary(f, rw, key, value);
        map->insert({key, value});
      }
    } else {
      for (auto &[key, value] : *map) {
        auto t_key = key;
        rw_binary(f, rw, t_key, value);
      }
    }
  };

  rw_map(S1.operator->());
  rw_map(dS1.operator->());
  rw_map(&e_sigma);
  rw_map(&en);

  // A file holding empty tables is treated as a failed read: downstream,
  // empty() means "no correction"
  return !readQ || !empty();
}

//==============================================================================
Sigma1Correction
calculate_dSdE_correction(const std::vector<DiracSpinor> &ci_basis,
                          const std::vector<DiracSpinor> &s1_basis_core,
                          const std::vector<DiracSpinor> &s1_basis_excited,
                          const Coulomb::QkTable &qk) {
  // nb: E0 is left at 0.0 (not set): it belongs to a single (J, parity), and
  // is found by iterate_E0()
  Sigma1Correction corr;

  for (const auto &v : ci_basis) {
    corr.en[v.nk_index()] = v.en();
    // First of each kappa: the energy Sigma_1 is evaluated at (same
    // convention as calculate_h1_table)
    corr.e_sigma.insert({v.kappa(), v.en()});
  }

  // First, in serial, add a zero for each element. The map structure is then
  // fixed, so the values may be filled in parallel below
  for (const auto &v : ci_basis) {
    for (const auto &w : ci_basis) {
      if (w.kappa() == v.kappa()) {
        corr.S1.add(v, w, 0.0);
        corr.dS1.add(v, w, 0.0);
      }
    }
  }

  qip::ProgressBar bar(ci_basis.size());
#pragma omp parallel for schedule(dynamic)
  for (auto iv = 0ul; iv < ci_basis.size(); ++iv) {
    const auto &v = ci_basis[iv];
    const auto ev = corr.e_sigma.at(v.kappa());

    for (const auto &w : ci_basis) {
      if (w > v)
        continue;
      if (w.kappa() != v.kappa())
        continue;

      const auto S0 =
        MBPT::Sigma_vw(v, w, qk, s1_basis_core, s1_basis_excited, 99, ev);
      const auto dS =
        MBPT::dSigma_dE_vw(v, w, qk, s1_basis_core, s1_basis_excited, ev);

      *corr.S1.get(v, w) = S0;
      *corr.dS1.get(v, w) = dS;
      // Add symmetric partner:
      if (v != w) {
        *corr.S1.get(w, v) = S0;
        *corr.dS1.get(w, v) = dS;
      }
    }
    bar.update();
  }
  return corr;
}

//==============================================================================
Sigma1Correction
calculate_dSdE_correction(const std::vector<DiracSpinor> &ci_basis,
                          const MBPT::CorrelationPotential &Sigma) {
  // nb: E0 is left at 0.0 (not set): it belongs to a single (J, parity), and
  // is found by iterate_E0()
  Sigma1Correction corr;

  for (const auto &v : ci_basis) {
    corr.en[v.nk_index()] = v.en();
    // e_sigma: the energy Sigma was actually formed at (first of each kappa;
    // same entry SigmaFv falls back to for basis states)
    const auto Sig = Sigma.get(v.kappa());
    if (Sig != nullptr) {
      corr.e_sigma.insert({v.kappa(), Sig->en});
    }
  }

  for (const auto &w : ci_basis) {
    // Identical to the Sigma_1 part of the h1 table (lambda + ladder
    // included); derivative matches the same Sigma entry (no ladder)
    const auto S_w = Sigma.SigmaFv(w);
    const auto dS_w = Sigma.dSigmaFv(w);

    for (const auto &v : ci_basis) {
      if (w > v)
        continue;
      if (w.kappa() != v.kappa())
        continue;

      const auto S0 = v * S_w;
      const auto dS = v * dS_w;
      corr.S1.add(v, w, S0);
      corr.dS1.add(v, w, dS);
      // Add symmetric partner:
      if (v != w) {
        corr.S1.add(w, v, S0);
        corr.dS1.add(w, v, dS);
      }
    }
  }
  return corr;
}

//==============================================================================
template <typename MatrixType>
MatrixType iterate_E0(PsiJPi *psi, const Sigma1Correction &s1c,
                      const Coulomb::meTable<double> &h1,
                      const Coulomb::QkTable &qk, const Coulomb::WkTable *Bk,
                      const Coulomb::WkTable *Sk, const std::vector<double> &hk,
                      const Coulomb::WkTable *dSk, double E0_sigma2,
                      std::ostream &outstream, std::size_t level,
                      const Coulomb::meTable<double> *Sigma1_in_CI) {
  assert(psi != nullptr);

  constexpr auto eps_E0 = 1.0e-6;
  constexpr auto max_iterations = 12;

  // Sigma_2 was tabulated at the reference E0_sigma2, and is shifted to the
  // current E0 using dS^k/dE0 (see corrected_Sk). Both corrections then move
  // together as E0 is iterated.
  const auto dE0_s2 = [dSk, E0_sigma2](double E0) {
    return dSk == nullptr ? 0.0 : E0 - E0_sigma2;
  };

  // E0 belongs to this (J, parity), but the correction tables are shared
  // between them: work on a local copy. Only E0 changes between passes.
  auto s1_Jpi = s1c;

  // E0 = 0.0 means "not set": start from the lowest zeroth-order configuration
  // (all valence electrons in the lowest orbital)
  if (s1_Jpi.E0 == 0.0) {
    auto e_min = 0.0;
    for (const auto &orbital : s1_Jpi.en) {
      if (orbital.second < e_min) {
        e_min = orbital.second;
      }
    }
    s1_Jpi.E0 = psi->num_electrons() * e_min;
  }

  fmt::print(outstream, "Iterating E0{} (initial: {:.6f}):\n",
             fmt::format(" ({})", level), s1_Jpi.E0);
  outstream << std::flush;

  // First pass: build Hci at the current E0, and diagonalise up to the
  // level. E0 enters only through the (small) dSigma/dE correction, so the
  // state itself hardly changes from pass to pass. Later passes therefore
  // rebuild Hci with the updated E0, and take the new E0 as the expectation
  // value of the new Hci in the (unchanged) state: first-order perturbation
  // theory in the change to Hci. Solving (dense: N^3; sparse: many products
  // with H) costs far more than one build, so this is much cheaper than
  // re-solving each pass.
  auto Hci = construct_Hci<MatrixType>(*psi, h1, qk, Bk, Sk, &s1_Jpi, hk, dSk,
                                       dE0_s2(s1_Jpi.E0), Sigma1_in_CI);
  if constexpr (std::is_same_v<MatrixType, LinAlg::Matrix<double>>) {
    psi->solve(Hci, int(level) + 1);
  } else {
    if (!psi->solve(Hci, int(level) + 1).converged) {
      fmt::print(outstream,
                 "Warning: Davidson did not converge (level {}, for E0)\n",
                 level);
    }
  }
  auto c_view = psi->coefs(level);
  const LinAlg::Vector<double> c(
    std::vector<double>(c_view.data(), c_view.data() + c_view.size()));

  auto delta = psi->energy(level) - s1_Jpi.E0;
  s1_Jpi.E0 = psi->energy(level);
  fmt::print(outstream, "  it {:2}: E0 = {:.6f} au (eps = {:.0e})\n", 1,
             s1_Jpi.E0, std::abs(delta / s1_Jpi.E0));
  outstream << std::flush;

  for (int it = 2; it <= max_iterations && std::abs(delta) >= eps_E0; ++it) {
    // Free the previous matrix first: only one is held at a time
    Hci = MatrixType{};
    Hci = construct_Hci<MatrixType>(*psi, h1, qk, Bk, Sk, &s1_Jpi, hk, dSk,
                                    dE0_s2(s1_Jpi.E0), Sigma1_in_CI);

    // <c|Hci|c>; c is normalised
    const double E0_new = c * (Hci * c);

    delta = E0_new - s1_Jpi.E0;
    s1_Jpi.E0 = E0_new;
    fmt::print(outstream, "  it {:2}: E0 = {:.8f} au (eps = {:.0e})\n", it,
               s1_Jpi.E0, std::abs(delta / s1_Jpi.E0));
    outstream << std::flush;
    if (it == max_iterations && std::abs(delta) >= eps_E0) {
      outstream << "Warning: E0 iteration did not converge\n";
    }
  }

  return Hci;
}

template LinAlg::Matrix<double> iterate_E0<LinAlg::Matrix<double>>(
  PsiJPi *, const Sigma1Correction &, const Coulomb::meTable<double> &,
  const Coulomb::QkTable &, const Coulomb::WkTable *, const Coulomb::WkTable *,
  const std::vector<double> &, const Coulomb::WkTable *, double, std::ostream &,
  std::size_t, const Coulomb::meTable<double> *);
template LinAlg::SparseSymmetricMatrix<double>
iterate_E0<LinAlg::SparseSymmetricMatrix<double>>(
  PsiJPi *, const Sigma1Correction &, const Coulomb::meTable<double> &,
  const Coulomb::QkTable &, const Coulomb::WkTable *, const Coulomb::WkTable *,
  const std::vector<double> &, const Coulomb::WkTable *, double, std::ostream &,
  std::size_t, const Coulomb::meTable<double> *);

//==============================================================================
// CI Hamiltonian matrix element for two 1-particle CSFs: one-body only. The
// Sigma_1 already in the CI space (hole CSFs in the block) is removed
double Hab(const CI::CSF1 &A, const CI::CSF1 &B, int,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &,
           const Sigma1Correction *s1c,
           const Coulomb::meTable<double> *Sigma1_in_CI) {
  const auto a = A.state(0);
  const auto b = B.state(0);
  return h1.getv(a, b) - (Sigma1_in_CI ? Sigma1_in_CI->getv(a, b) : 0.0) +
         (s1c ? s1c->delta_h1(a, b) : 0.0);
}

//==============================================================================
double CSF2_one_body(const Coulomb::meTable<double> &h1, DiracSpinor::Index v,
                     DiracSpinor::Index w, DiracSpinor::Index x,
                     DiracSpinor::Index y, int twoJ,
                     const Sigma1Correction *s1c) {

  const auto etaV = v == w ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto etaX = x == y ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto eta2 = etaV * etaX;

  const auto tjv = Angular::nkindex_to_twoj(v);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto s = Angular::neg1pow_2(twoJ + tjv + tjw);

  // In each term, the matched orbital is the spectator (its energy enters
  // the dSigma/dE correction, if included)
  double h1_VX = 0.0;
  if (x == v) {
    h1_VX += eta2 * (h1.getv(y, w) + (s1c ? s1c->delta_h1(y, w, v) : 0.0));
  }
  if (x == w) {
    h1_VX -= s * eta2 * (h1.getv(y, v) + (s1c ? s1c->delta_h1(y, v, w) : 0.0));
  }
  if (y == v) {
    h1_VX -= s * eta2 * (h1.getv(x, w) + (s1c ? s1c->delta_h1(x, w, v) : 0.0));
  }
  if (y == w) {
    h1_VX += eta2 * (h1.getv(x, v) + (s1c ? s1c->delta_h1(x, v, w) : 0.0));
  }
  return h1_VX;
}

//==============================================================================
// Determines CI Hamiltonian matrix element for two 2-particle CSFs, a and b
double Hab(const CI::CSF2 &X, const CI::CSF2 &V, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
           const Sigma1Correction *s1c) {
  const auto [v, w] = V.states;
  const auto [x, y] = X.states;
  return CSF2_one_body(h1, v, w, x, y, twoJ, s1c) +
         CSF2_Coulomb(qk, v, w, x, y, twoJ);
}

//==============================================================================
PairElement coulomb_pair_element(const Coulomb::QkTable &qk) {
  return [&qk](DiracSpinor::Index v, DiracSpinor::Index w, DiracSpinor::Index x,
               DiracSpinor::Index y,
               int twoJ) { return CSF2_Coulomb(qk, v, w, x, y, twoJ); };
}

//==============================================================================
PairElement breit_pair_element(const Coulomb::WkTable &Bk) {
  return [&Bk](DiracSpinor::Index v, DiracSpinor::Index w, DiracSpinor::Index x,
               DiracSpinor::Index y,
               int twoJ) { return CSF2_Breit(Bk, v, w, x, y, twoJ); };
}

//==============================================================================
PairElement sigma2_pair_element(const Coulomb::WkTable &Sk,
                                const Coulomb::QkTable *qk,
                                const std::vector<double> &hk,
                                const Coulomb::WkTable *dSk, double dE0) {
  return [&Sk, qk, &hk, dSk, dE0](DiracSpinor::Index v, DiracSpinor::Index w,
                                  DiracSpinor::Index x, DiracSpinor::Index y,
                                  int twoJ) {
    return CSF2_Sigma2(Sk, v, w, x, y, twoJ, qk, hk, dSk, dE0);
  };
}

//==============================================================================
double particle_hole_element(const PairElement &pair_element,
                             DiracSpinor::Index y, DiracSpinor::Index b,
                             DiracSpinor::Index w, DiracSpinor::Index a,
                             int twoK) {
  const auto tjy = Angular::nkindex_to_twoj(y);
  const auto tjb = Angular::nkindex_to_twoj(b);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tja = Angular::nkindex_to_twoj(a);

  // Recoupling: the element between the particle-particle pairs (y a) and
  // (w b), coupled to J', summed with a 6j symbol over J' to give the element
  // between the particle-hole pairs (y b~) and (w a~) coupled to K.
  // The overall minus sign comes from reordering the hole operators
  const auto twoJ_min = std::max(std::abs(tjy - tja), std::abs(tjw - tjb));
  const auto twoJ_max = std::min(tjy + tja, tjw + tjb);
  double out = 0.0;
  for (int twoJp = twoJ_min; twoJp <= twoJ_max; twoJp += 2) {
    const auto sj = Angular::sixj_2(tjw, tja, twoK, tjy, tjb, twoJp);
    if (sj == 0.0)
      continue;
    out -= (twoJp + 1) * sj * pair_element(y, a, w, b, twoJp);
  }
  return out;
}

//==============================================================================
double two_body_2p1h(const PairElement &pair_element, const CI::CSF2p1h &A,
                     const CI::CSF2p1h &B, int twoJ) {

  if (CSF2p1h::num_different(A, B) > 2)
    return 0.0;

  const auto [x, y] = A.states;
  const auto [v, w] = B.states;

  double out = 0.0;

  // Between the particles: hole spectator
  if (A.hole == B.hole && A.twoJ12() == B.twoJ12()) {
    out += pair_element(v, w, x, y, B.twoJ12());
  }

  // Between each particle and the hole
  const auto ph_element =
    [&pair_element](DiracSpinor::Index yp, DiracSpinor::Index bp,
                    DiracSpinor::Index wp, DiracSpinor::Index ap, int twoK) {
      return particle_hole_element(pair_element, yp, bp, wp, ap, twoK);
    };
  out += particle_hole_sum(A, B, twoJ, ph_element);

  return out;
}

//==============================================================================
double two_body_1p_2p1h(const PairElement &pair_element, const CI::CSF1 &A,
                        const CI::CSF2p1h &B, int twoJ) {
  const auto x = A.state(0);
  const auto [v, w] = B.states;
  const auto a = B.hole;
  const auto tjx = Angular::nkindex_to_twoj(x);
  const auto tja = Angular::nkindex_to_twoj(a);
  const auto twoJ12 = B.twoJ12();
  assert(tjx == twoJ);
  // The pair (vw), coupled to J12, de-excites into the electron x and the
  // hole: the pair element <(xa)J12|g~|(vw)J12> (which includes eta_vw)
  return Angular::neg1pow_2(tjx - tja + twoJ12) *
         std::sqrt(double(twoJ12 + 1) / (twoJ + 1)) *
         pair_element(x, a, v, w, twoJ12);
}

//==============================================================================
double Hab(const CI::CSF2p1h &A, const CI::CSF2p1h &B, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
           const Sigma1Correction *) {

  if (CSF2p1h::num_different(A, B) > 2)
    return 0.0;

  const auto [x, y] = A.states;
  const auto [v, w] = B.states;
  const auto b = A.hole;
  const auto a = B.hole;
  const auto same_J12 = A.twoJ12() == B.twoJ12();

  double out = 0.0;

  // One-body on the particles: hole spectator
  if (a == b && same_J12) {
    out += CSF2_one_body(h1, v, w, x, y, B.twoJ12());
  }

  // One-body on the hole: pair spectator
  if (A.states == B.states && same_J12) {
    out -= h1.getv(a, b);
  }

  // Coulomb between the particles, and between each particle and the hole
  out += two_body_2p1h(coulomb_pair_element(qk), A, B, twoJ);

  return out;
}

//==============================================================================
double Hab(const CI::CSF1 &A, const CI::CSF2p1h &B, int twoJ,
           const Coulomb::meTable<double> &, const Coulomb::QkTable &qk,
           const Sigma1Correction *) {
  // The one-body part vanishes: <w|h_HF|a> = 0 for a core orbital a
  return two_body_1p_2p1h(coulomb_pair_element(qk), A, B, twoJ);
}

//==============================================================================
double Hab(const CI::CSF3 &A, const CI::CSF3 &B, int,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
           const Sigma1Correction *s1c) {

  if (CSF3::num_different(A, B) > 2)
    return 0.0;

  double out = 0.0;
  for (const auto &p : A.parents()) {
    for (const auto &q : B.parents()) {
      if (p.twoJ2 != q.twoJ2)
        continue;
      const auto c = p.coef * q.coef;

      // One-body: acts on electron 3; the parent pair must coincide, and its
      // two orbitals are the spectators for the dSigma/dE correction
      if (p.pair == q.pair) {
        const auto [x, y] = p.pair;
        const auto dh1 =
          s1c ? s1c->delta_h1(p.spectator, q.spectator, x, y) : 0.0;
        out += c * (h1.getv(p.spectator, q.spectator) + dh1);
      }

      // Two-body: acts on electrons 1 and 2; the third electron must coincide
      if (p.spectator == q.spectator) {
        out += c * CSF2_Coulomb(qk, p.pair[0], p.pair[1], q.pair[0], q.pair[1],
                                p.twoJ2);
      }
    }
  }
  return 3.0 * out;
}

//==============================================================================
double Hab(const CI::CSF &A, const CI::CSF &B, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
           const Sigma1Correction *s1c,
           const Coulomb::meTable<double> *Sigma1_in_CI) {
  return std::visit(
    [&](const auto &a, const auto &b) -> double {
      using A_type = std::decay_t<decltype(a)>;
      using B_type = std::decay_t<decltype(b)>;
      if constexpr (std::is_same_v<A_type, CSF1> &&
                    std::is_same_v<B_type, CSF1>) {
        return Hab(a, b, twoJ, h1, qk, s1c, Sigma1_in_CI);
      } else if constexpr (std::is_same_v<A_type, B_type>) {
        return Hab(a, b, twoJ, h1, qk, s1c);
      } else if constexpr (std::is_same_v<A_type, CSF1> &&
                           std::is_same_v<B_type, CSF2p1h>) {
        return Hab(a, b, twoJ, h1, qk, s1c);
      } else if constexpr (std::is_same_v<A_type, CSF2p1h> &&
                           std::is_same_v<B_type, CSF1>) {
        // Real and symmetric
        return Hab(b, a, twoJ, h1, qk, s1c);
      } else {
        assert(false && "CI matrix elements between CSFs of different kinds "
                        "are not implemented");
        return 0.0;
      }
    },
    A.variant(), B.variant());
}

//==============================================================================
Coulomb::meTable<double>
calculate_h1_table(const std::vector<DiracSpinor> &ci_basis,
                   const std::vector<DiracSpinor> &s1_basis_core,
                   const std::vector<DiracSpinor> &s1_basis_excited,
                   const Coulomb::QkTable &qk, bool include_Sigma1) {
  // Create lookup table for one-particle matrix elements, h1
  Coulomb::meTable<double> h1;

  // First, in serial, add a zero for each element. The map structure is then
  // fixed, so the values may be filled in parallel below
  for (const auto &v : ci_basis) {
    for (const auto &w : ci_basis) {
      if (w.kappa() == v.kappa())
        h1.add(v, w, 0.0);
    }
  }

// Calculate + store all 1-body integrals
// Each v writes only its own elements, so this is safe in parallel
#pragma omp parallel for schedule(dynamic)
  for (auto iv = 0ul; iv < ci_basis.size(); ++iv) {
    const auto &v = ci_basis[iv];

    // Find lowest valence state of this kappa; for Sigma energy
    const auto vp = std::find_if(
      ci_basis.cbegin(), ci_basis.cend(),
      [kappa = v.kappa()](const auto &n) { return n.kappa() == kappa; });
    auto ev = vp->en();

    for (const auto &w : ci_basis) {
      if (w > v)
        continue;
      if (w.kappa() != v.kappa())
        continue;
      const auto h0_vw = v == w ? v.en() : 0.0;

      const auto Sigma_vw =
        include_Sigma1 ?
          MBPT::Sigma_vw(v, w, qk, s1_basis_core, s1_basis_excited, 99, ev) :
          0.0;

      *h1.get(v, w) = h0_vw + Sigma_vw;
      // Add symmetric partner:
      if (v != w)
        *h1.get(w, v) = h0_vw + Sigma_vw;
    }
  }
  return h1;
}

//==============================================================================
Coulomb::meTable<double>
calculate_Sigma1_in_CI(const std::vector<DiracSpinor> &ci_basis,
                       const std::vector<DiracSpinor> &holes,
                       const std::vector<DiracSpinor> &holes_basis,
                       const ConfigurationSpace &config_space,
                       const Coulomb::QkTable &qk) {
  Coulomb::meTable<double> Sigma1_in_CI;

  // The intermediate states (b, m, n) that are hole CSFs of the block
  const auto in_CI_space = [&config_space](const DiracSpinor &b,
                                           const DiracSpinor &m,
                                           const DiracSpinor &n) {
    return config_space.included({m.nk_index(), n.nk_index()}, {b.nk_index()});
  };

  for (const auto &v : ci_basis) {
    // Same energy as calculate_h1_table: the lowest orbital of this kappa
    const auto vp = std::find_if(
      ci_basis.cbegin(), ci_basis.cend(),
      [kappa = v.kappa()](const auto &n) { return n.kappa() == kappa; });
    const auto ev = vp->en();

    for (const auto &w : ci_basis) {
      if (w > v || w.kappa() != v.kappa())
        continue;
      const auto Sigma_P =
        MBPT::Sigma_vw_2p1h(v, w, qk, holes, holes_basis, ev, in_CI_space);
      Sigma1_in_CI.add(v, w, Sigma_P);
      if (v != w)
        Sigma1_in_CI.add(w, v, Sigma_P);
    }
  }
  return Sigma1_in_CI;
}

//==============================================================================
Coulomb::meTable<double>
calculate_h1_table(const std::vector<DiracSpinor> &ci_basis,
                   const MBPT::CorrelationPotential &Sigma,
                   bool include_Sigma1) {
  // Create lookup table for one-particle matrix elements, h1
  Coulomb::meTable<double> h1;

  // Calculate + store all 1-body integrals
  for (const auto &w : ci_basis) {

    const auto Sigma_w = Sigma(w);

    for (const auto &v : ci_basis) {
      if (w > v)
        continue;
      if (w.kappa() != v.kappa())
        continue;
      const auto h0_vw = v == w ? v.en() : 0.0;

      // Can use Sigma matrix instead: all-orders?
      const auto Sigma_vw = include_Sigma1 ? v * Sigma_w : 0.0;

      h1.add(v, w, h0_vw + Sigma_vw);
      // Add symmetric partner:
      if (v != w)
        h1.add(w, v, h0_vw + Sigma_vw);
    }
  }
  return h1;
}

//==============================================================================
Coulomb::WkTable calculate_Bk(const std::string &bk_filename,
                              const HF::Breit *const pBr,
                              const std::vector<DiracSpinor> &ci_basis,
                              int max_k, bool no_new_integrals) {
  // Breit table
  Coulomb::WkTable Bk;

  Bk.read(bk_filename);
  const auto existing = Bk.count();

  if (!no_new_integrals) {

    auto vBr = *pBr;              // copy, so we can fill two-particle 'bk'
    vBr.fill_gb(ci_basis, max_k); // very quick, and makes below *much* faster
    //* nb: uses HUGE memory if basis is large; CI basis normally small enough

    const auto Bk_function = [&](int k, const DiracSpinor &v,
                                 const DiracSpinor &w, const DiracSpinor &x,
                                 const DiracSpinor &y) {
      return vBr.Bk_abcd_2(k, v, w, x, y);
    };

    Bk.fill(ci_basis, Bk_function, HF::Breit::Bk_SR, max_k, false);

    // print summary
    Bk.summary();

    // If we calculated new integrals, write to disk
    const auto total = Bk.count();
    assert(total >= existing);
    const auto new_integrals = total - existing;
    std::cout << "Calculated " << new_integrals << " new Breit integrals\n";
    if (new_integrals > 0) {
      Bk.write(bk_filename);
    }
    std::cout << "\n" << std::flush;
  }

  return Bk;
}

//==============================================================================
// Takes a subset of input basis according to subset_string.
// Only states *not* included in frozen_core_string are included.
std::vector<DiracSpinor> basis_subset(const std::vector<DiracSpinor> &basis,
                                      const std::string &subset_string,
                                      const std::string &frozen_core_string) {

  // Form 'subset' from {a} in 'basis', if:
  //    a _is_ in subset_string (or subset_string is empty) AND
  //    a _is not_ in frozen_core_string

  std::vector<DiracSpinor> subset;
  const auto nmaxk_list = AtomData::n_kappa_list(subset_string);
  const auto core_list = AtomData::core_parser(frozen_core_string);

  for (const auto &a : basis) {

    // If subset_string is non-empty, check that a is within it
    if (!subset_string.empty()) {
      const auto nk =
        std::find_if(nmaxk_list.cbegin(), nmaxk_list.cend(),
                     [&a](const auto &tnk) { return a.kappa() == tnk.second; });
      if (nk == nmaxk_list.cend())
        continue;
      // n labels of negative-energy states continue on from the electron
      // states, so the n limit applies to the electron states only
      if (!a.negativeEnergyStateQ() && a.n() > nk->first)
        continue;
    }

    // assume only filled shells in frozen core
    const auto core = std::find_if(
      core_list.cbegin(), core_list.cend(),
      [&a](const auto &tcore) { return a.n() == tcore.n && a.l() == tcore.l; });

    if (core != core_list.cend())
      continue;
    subset.push_back(a);
  }
  return subset;
}

//==============================================================================
// Calculate reduced matrix elements between two CI states.
// cA is CI expansion coefficients (row if CI eigenvector matrix)
double ReducedME(const LinAlg::View<const double> &cA,
                 const std::vector<CSF> &CSFAs, int twoJA,
                 const LinAlg::View<const double> &cB,
                 const std::vector<CSF> &CSFBs, int twoJB,
                 const Coulomb::meTable<double> &h, int K_rank, int Parity) {

  // selection rules: Not required (operator encodes it's own),
  // but it's faster to check in advance
  const auto piA = CSFAs.front().parity();
  const auto piB = CSFBs.front().parity();
  if (piA * piB != Parity)
    return 0.0;
  if (std::abs(twoJA - twoJB) > 2 * K_rank)
    return 0.0;
  if (K_rank != 0 && (twoJA == 0 && twoJB == 0))
    return 0.0;

  // <A|h|A>   = Σ_{IJ} c_I * c_J * <I|h|J>
  // <A||h||A> = Σ_{IJ} c_I * c_J * <I||h||J>
  // Only pairs of CSFs differing by at most one orbital contribute
  const auto connections = one_body_connected(CSFAs, CSFBs);

  double rme = 0.0;
#pragma omp parallel for schedule(dynamic) reduction(+ : rme)
  for (std::size_t a = 0ul; a < connections.size(); ++a) {
    for (const auto i : connections[a].config_CSFs) {
      const auto &csf_i = CSFAs[i];
      const auto ci = cA[i];
      for (const auto j : connections[a].connected_CSFs) {
        const auto &csf_j = CSFBs[j];
        const auto cj = cB[j];

        rme += ci * cj * RME_CSF(csf_i, twoJA, csf_j, twoJB, h, K_rank);
      }
    }
  }
  return rme;
}

//==============================================================================
double RME_CSF1(const CI::CSF1 &X, int, const CI::CSF1 &V, int,
                const Coulomb::meTable<double> &h, int) {
  return h.getv(X.state(0), V.state(0));
}

//==============================================================================
double RME_CSF3(const CI::CSF3 &X, int twoJX, const CI::CSF3 &V, int twoJV,
                const Coulomb::meTable<double> &h, int K_rank) {

  if (CSF3::num_different(X, V) > 1)
    return 0.0;

  const auto twok = 2 * K_rank;
  const auto f = 3.0 * std::sqrt(double(twoJX + 1) * (twoJV + 1));

  double sum = 0.0;
  for (const auto &p : X.parents()) {
    for (const auto &q : V.parents()) {
      // One-body: acts on electron 3; the parent pair must coincide
      if (p.twoJ2 != q.twoJ2 || p.pair != q.pair)
        continue;
      const auto t = h.getv(p.spectator, q.spectator);
      if (t == 0.0)
        continue;
      const auto tjp = Angular::nkindex_to_twoj(p.spectator);
      const auto tjq = Angular::nkindex_to_twoj(q.spectator);
      const auto sj = Angular::sixj_2(tjp, twoJX, p.twoJ2, twoJV, tjq, twok);
      const auto s = Angular::neg1pow_2(p.twoJ2 + tjq + twoJX + twok);
      sum += p.coef * q.coef * s * sj * t;
    }
  }
  return f * sum;
}

//==============================================================================
// Calculate reduce ME between two 2-particle CSFs
double RME_CSF2(const CI::CSF2 &X, int twoJX, const CI::CSF2 &V, int twoJV,
                const Coulomb::meTable<double> &h, int K_rank) {
  return RME_pair(X.state(0), X.state(1), twoJX, V.state(0), V.state(1), twoJV,
                  h, K_rank);
}

//==============================================================================
double RME_pair(DiracSpinor::Index x, DiracSpinor::Index y, int twoJX,
                DiracSpinor::Index v, DiracSpinor::Index w, int twoJV,
                const Coulomb::meTable<double> &h, int K_rank) {

  const auto etaV = v == w ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto etaX = x == y ? 1.0 / std::sqrt(2.0) : 1.0;

  const auto twok = 2 * K_rank;

  const auto f = etaV * etaX * std::sqrt(double(twoJX + 1) * (twoJV + 1)) *
                 Angular::neg1pow(K_rank);

  const auto tjv = Angular::nkindex_to_twoj(v);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tjx = Angular::nkindex_to_twoj(x);
  const auto tjy = Angular::nkindex_to_twoj(y);

  double sum = 0.0;
  if (y == w) {
    const auto sj = Angular::sixj_2(twoJX, twoJV, twok, tjv, tjx, tjw);
    const auto t = h.getv(x, v);
    const auto s = Angular::neg1pow_2(tjw + tjx + twoJV);
    sum += f * sj * t * s;
  }
  if (y == v) {
    const auto sj = Angular::sixj_2(twoJX, twoJV, twok, tjw, tjx, tjv);
    const auto t = h.getv(x, w);
    const auto s = Angular::neg1pow_2(tjw + tjx);
    sum += f * sj * t * s;
  }
  if (x == w) {
    const auto sj = Angular::sixj_2(twoJX, twoJV, twok, tjv, tjy, tjw);
    const auto t = h.getv(y, v);
    const auto s = Angular::neg1pow_2(twoJX + twoJV + 2);
    sum += f * sj * t * s;
  }
  if (x == v) {
    const auto sj = Angular::sixj_2(twoJX, twoJV, twok, tjw, tjy, tjv);
    const auto t = h.getv(y, w);
    const auto s = Angular::neg1pow_2(tjw + tjv + twoJX);
    sum += f * sj * t * s;
  }
  return sum;
}

//==============================================================================
double RME_CSF2p1h(const CI::CSF2p1h &X, int twoJX, const CI::CSF2p1h &V,
                   int twoJV, const Coulomb::meTable<double> &h, int K_rank) {

  if (CSF2p1h::num_different(X, V) > 1)
    return 0.0;

  const auto twok = 2 * K_rank;
  const auto b = X.hole;
  const auto a = V.hole;
  const auto tja = Angular::nkindex_to_twoj(a);
  const auto tjb = Angular::nkindex_to_twoj(b);
  const auto twoJ12_X = X.twoJ12();
  const auto twoJ12_V = V.twoJ12();
  const auto f = std::sqrt(double(twoJX + 1) * (twoJV + 1));

  double out = 0.0;

  // On a particle: the pair is the first part of the coupling, the hole is
  // the spectator
  if (a == b) {
    const auto sj =
      Angular::sixj_2(twoJ12_X, twoJX, tja, twoJV, twoJ12_V, twok);
    if (sj != 0.0) {
      const auto s = Angular::neg1pow_2(twoJ12_X + tja + twoJV + twok);
      out += s * f * sj *
             RME_pair(X.state(0), X.state(1), twoJ12_X, V.state(0), V.state(1),
                      twoJ12_V, h, K_rank);
    }
  }

  // On the hole: the pair is the spectator. Hole single-particle element:
  // <b~||t||a~> = (-1)^(ja + jb + k) <a||t||b>
  if (X.states == V.states && twoJ12_X == twoJ12_V) {
    const auto t_ab = h.getv(a, b);
    const auto sj = Angular::sixj_2(tjb, twoJX, twoJ12_V, twoJV, tja, twok);
    if (t_ab != 0.0 && sj != 0.0) {
      const auto s = Angular::neg1pow_2(twoJ12_V + tja + twoJX + twok) *
                     Angular::neg1pow_2(tja + tjb + twok);
      out += s * f * sj * t_ab;
    }
  }
  return out;
}

//==============================================================================
double RME_CSF2p1h_1(const CI::CSF2p1h &X, int twoJX, const CI::CSF1 &V,
                     int twoJV, const Coulomb::meTable<double> &h, int K_rank) {
  const auto x = V.state(0);
  const auto [v, w] = X.states;
  const auto a = X.hole;
  if (x != v && x != w)
    return 0.0;

  const auto twok = 2 * K_rank;
  const auto tjv = Angular::nkindex_to_twoj(v);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tja = Angular::nkindex_to_twoj(a);
  assert(Angular::nkindex_to_twoj(x) == twoJV);
  const auto twoJ12 = X.twoJ12();
  const auto eta = v == w ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto f = -Angular::neg1pow_2(tjv + tjw + tja + twoJX) *
                 std::sqrt(double(twoJX + 1) * (twoJ12 + 1)) * eta;
  const auto s_vw = Angular::neg1pow_2(tjv + tjw - twoJ12 + 2);

  // The electron x is one of the two particles; the operator excites the
  // hole into the other
  double out = 0.0;
  if (x == v) {
    out += Angular::sixj_2(tjv, tjw, twoJ12, tja, twoJX, twok) * h.getv(w, a);
  }
  if (x == w) {
    out +=
      s_vw * Angular::sixj_2(tjw, tjv, twoJ12, tja, twoJX, twok) * h.getv(v, a);
  }
  return f * out;
}

//==============================================================================
double RME_CSF1_2p1h(const CI::CSF1 &X, int twoJX, const CI::CSF2p1h &V,
                     int twoJV, const Coulomb::meTable<double> &h, int K_rank) {
  const auto x = X.state(0);
  const auto [v, w] = V.states;
  const auto a = V.hole;
  if (x != v && x != w)
    return 0.0;

  const auto twok = 2 * K_rank;
  const auto tjv = Angular::nkindex_to_twoj(v);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tja = Angular::nkindex_to_twoj(a);
  const auto tjx = Angular::nkindex_to_twoj(x);
  assert(tjx == twoJX);
  const auto twoJ12 = V.twoJ12();
  const auto eta = v == w ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto f = -Angular::neg1pow_2(tjv + tjw + tja + twoJV) *
                 std::sqrt(double(twoJV + 1) * (twoJ12 + 1)) * eta;
  const auto s_vw = Angular::neg1pow_2(tjv + tjw - twoJ12 + 2);

  // The electron x is one of the two particles; the operator de-excites the
  // other into the hole (the conjugate of RME_CSF2p1h_1, term by term)
  double out = 0.0;
  if (x == v) {
    out += Angular::neg1pow_2(tja - tjw + twoJV - tjx) *
           Angular::sixj_2(tjv, tjw, twoJ12, tja, twoJV, twok) * h.getv(a, w);
  }
  if (x == w) {
    out += s_vw * Angular::neg1pow_2(tja - tjv + twoJV - tjx) *
           Angular::sixj_2(tjw, tjv, twoJ12, tja, twoJV, twok) * h.getv(a, v);
  }
  return f * out;
}

//==============================================================================
double RME_CSF(const CI::CSF &X, int twoJX, const CI::CSF &V, int twoJV,
               const Coulomb::meTable<double> &h, int K_rank) {
  // Reduced matrix element for the kinds of CSF held by X and V
  const auto rme_by_kind = [&](const auto &x, const auto &v) -> double {
    using X_type = std::decay_t<decltype(x)>;
    using V_type = std::decay_t<decltype(v)>;
    if constexpr (std::is_same_v<X_type, CSF1> &&
                  std::is_same_v<V_type, CSF2p1h>) {
      return RME_CSF1_2p1h(x, twoJX, v, twoJV, h, K_rank);
    } else if constexpr (std::is_same_v<X_type, CSF2p1h> &&
                         std::is_same_v<V_type, CSF1>) {
      return RME_CSF2p1h_1(x, twoJX, v, twoJV, h, K_rank);
    } else if constexpr (!std::is_same_v<X_type, V_type>) {
      assert(false && "Reduced matrix elements between CSFs of different "
                      "kinds are not implemented");
      return 0.0;
    } else if constexpr (std::is_same_v<X_type, CSF1>) {
      return RME_CSF1(x, twoJX, v, twoJV, h, K_rank);
    } else if constexpr (std::is_same_v<X_type, CSF2>) {
      return RME_CSF2(x, twoJX, v, twoJV, h, K_rank);
    } else if constexpr (std::is_same_v<X_type, CSF2p1h>) {
      return RME_CSF2p1h(x, twoJX, v, twoJV, h, K_rank);
    } else {
      return RME_CSF3(x, twoJX, v, twoJV, h, K_rank);
    }
  };

  return std::visit(rme_by_kind, X.variant(), V.variant());
}

//==============================================================================
std::pair<std::string, double>
leading_config(const LinAlg::View<const double> &coefs,
               const std::vector<CSF> &csfs) {
  assert(coefs.size() == csfs.size());

  // Total |c|^2 of each non-rel configuration (summed over its rel. CSFs):
  std::map<std::string, double> nr_configs{};
  for (std::size_t j = 0; j < csfs.size(); ++j) {
    nr_configs[csfs.at(j).config(false)] += coefs[j] * coefs[j];
  }

  const auto max_it = std::max_element(
    nr_configs.begin(), nr_configs.end(),
    [](const auto &a, const auto &b) { return a.second < b.second; });

  return max_it != nr_configs.end() ? std::pair{max_it->first, max_it->second} :
                                      std::pair{std::string{}, 0.0};
}

//==============================================================================
std::pair<int, int> Term_S_L(int l1, int l2, int twoJ, double gJ_target) {
  // Determine Term Symbol, from g-factor
  const auto min_L = std::abs(l1 - l2);
  const auto max_L = std::abs(l1 + l2);
  const auto min_S = 0;
  const auto max_S = 1;

  int L = -1;
  int S = min_S;
  double best_del = 2.0;

  for (int tL = min_L; tL <= max_L; ++tL) {
    for (int tS = min_S; tS <= max_S; ++tS) {

      if ((twoJ > 2 * (tL + tS)) || (twoJ < 2 * (tL - tS)))
        continue;

      if (twoJ == 0) {
        if (tL > L) {
          // choose smallest S for Largest allowed L
          L = tL;
          S = tS;
        }
      }

      if (twoJ != 0) {
        const auto gJNR = 1.5 + (tS * (tS + 1.0) - tL * (tL + 1.0)) /
                                  (twoJ * (0.5 * twoJ + 1.0));
        if (std::abs(gJ_target - gJNR) < best_del) {
          best_del = std::abs(gJ_target - gJNR);
          L = tL;
          S = tS;
        }
      }
    }
  }

  return {S, L};
}

//==============================================================================
std::pair<int, int> Term_S_L_from_expectation(double L2, double S2, int twoJ) {
  // Nearest (S, L) consistent with J. S has the parity of J: for two
  // electrons S = 0, 1; for three, S = 1/2, 3/2.
  // L limit: generous; sum of l's for any realistic CI basis
  constexpr int max_L = 16;

  int best_L = -1;
  int best_twoS = twoJ % 2;
  double best_del = 1.0e30;
  for (int twoS = twoJ % 2; twoS <= twoJ % 2 + 2; twoS += 2) {
    const auto S = 0.5 * twoS;
    for (int L = 0; L <= max_L; ++L) {
      if (twoJ > 2 * L + twoS || twoJ < std::abs(2 * L - twoS))
        continue;
      const auto del =
        std::abs(L * (L + 1.0) - L2) + std::abs(S * (S + 1.0) - S2);
      if (del < best_del) {
        best_del = del;
        best_L = L;
        best_twoS = twoS;
      }
    }
  }
  return {best_twoS, best_L};
}

//==============================================================================
std::string Term_Symbol(int two_J, int L, int two_S, int parity) {
  return two_J % 2 == 0 ?
           fmt::format("{}^{}{}_{}", two_S + 1, AtomData::L_symbol(L),
                       parity == 1 ? "" : "°", two_J / 2) :
           fmt::format("{}^{}{}_{}/2", two_S + 1, AtomData::L_symbol(L),
                       parity == 1 ? "" : "°", two_J);
}

std::string Term_Symbol(int L, int two_S, int parity) {
  return fmt::format("{}{}{}", two_S + 1, AtomData::L_symbol(L),
                     parity == 1 ? "" : "°");
}

//==============================================================================
template <typename MatrixType>
MatrixType
construct_Hci(const PsiJPi &psi, const Coulomb::meTable<double> &h1,
              const Coulomb::QkTable &qk, const Coulomb::WkTable *Bk,
              const Coulomb::WkTable *Sk, const Sigma1Correction *s1c,
              const std::vector<double> &hk, const Coulomb::WkTable *dSk,
              double dE0, const Coulomb::meTable<double> *Sigma1_in_CI) {
  return construct_Hci<MatrixType>(psi.CSFs(), psi.twoJ(), h1, qk, Bk, Sk, s1c,
                                   hk, dSk, dE0, Sigma1_in_CI);
}

template LinAlg::Matrix<double> construct_Hci<LinAlg::Matrix<double>>(
  const PsiJPi &, const Coulomb::meTable<double> &, const Coulomb::QkTable &,
  const Coulomb::WkTable *, const Coulomb::WkTable *, const Sigma1Correction *,
  const std::vector<double> &, const Coulomb::WkTable *, double,
  const Coulomb::meTable<double> *);
template LinAlg::SparseSymmetricMatrix<double>
construct_Hci<LinAlg::SparseSymmetricMatrix<double>>(
  const PsiJPi &, const Coulomb::meTable<double> &, const Coulomb::QkTable &,
  const Coulomb::WkTable *, const Coulomb::WkTable *, const Sigma1Correction *,
  const std::vector<double> &, const Coulomb::WkTable *, double,
  const Coulomb::meTable<double> *);

//==============================================================================
template <typename MatrixType>
MatrixType
construct_Hci(const std::vector<CSF> &CSFs, int twoJ,
              const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
              const Coulomb::WkTable *Bk, const Coulomb::WkTable *Sk,
              const Sigma1Correction *s1c, const std::vector<double> &hk,
              const Coulomb::WkTable *dSk, double dE0,
              const Coulomb::meTable<double> *Sigma1_in_CI) {

  constexpr bool sparse =
    std::is_same_v<MatrixType, LinAlg::SparseSymmetricMatrix<double>>;

  const auto N_CSFs = CSFs.size();
  MatrixType Hci(N_CSFs);

  // Row iA has iA + 1 elements: dynamic schedule balances the triangle
#pragma omp parallel for schedule(dynamic)
  for (std::size_t iA = 0; iA < N_CSFs; ++iA) {
    // Sparse: the non-zero elements of this row (lower triangle)
    std::vector<std::uint32_t> columns;
    std::vector<double> values;
    // go to iB <= iA only: symmetric matrix
    for (std::size_t iB = 0; iB <= iA; ++iB) {

      const auto &A = CSFs[iA];
      const auto &B = CSFs[iB];

      // CSFs differing by more than two orbitals do not interact (three
      // electrons: most pairs); skip before any element is evaluated
      if (num_different(A, B) > 2)
        continue;

      // Regular CI matrix (h1 may include Sigma_1):
      const auto E_AB = Hab(A, B, twoJ, h1, qk, s1c, Sigma1_in_CI);
      // Two-body corrections: Sigma_2 (three electrons: the pair energies
      // need the orbital energy of the third electron), and Breit
      const auto dEs_AB = Sk != nullptr ?
                            CI::Sigma2_AB(A, B, twoJ, *Sk, &qk, hk, dSk, dE0,
                                          s1c ? &s1c->en : nullptr) :
                            0.0;
      const auto dEb_AB = Bk != nullptr ? CI::Breit_AB(A, B, twoJ, *Bk) : 0.0;

      const auto H_AB = E_AB + dEs_AB + dEb_AB;
      if constexpr (sparse) {
        // Exact zeros only (e.g., CSFs differing by more than two orbitals)
        if (H_AB != 0.0) {
          columns.push_back(std::uint32_t(iB));
          values.push_back(H_AB);
        }
      } else {
        Hci(iA, iB) = H_AB;
        // fill other half of symmetric matrix:
        if (iB != iA) {
          Hci(iB, iA) = H_AB;
        }
      }
    }
    if constexpr (sparse) {
      Hci.set_row(iA, columns, values);
    }
  }

  return Hci;
}

template LinAlg::Matrix<double> construct_Hci<LinAlg::Matrix<double>>(
  const std::vector<CSF> &, int, const Coulomb::meTable<double> &,
  const Coulomb::QkTable &, const Coulomb::WkTable *, const Coulomb::WkTable *,
  const Sigma1Correction *, const std::vector<double> &,
  const Coulomb::WkTable *, double, const Coulomb::meTable<double> *);
template LinAlg::SparseSymmetricMatrix<double>
construct_Hci<LinAlg::SparseSymmetricMatrix<double>>(
  const std::vector<CSF> &, int, const Coulomb::meTable<double> &,
  const Coulomb::QkTable &, const Coulomb::WkTable *, const Coulomb::WkTable *,
  const Sigma1Correction *, const std::vector<double> &,
  const Coulomb::WkTable *, double, const Coulomb::meTable<double> *);

//==============================================================================
LinAlg::Matrix<double> construct_Hci(const PsiJPi &psi, const Integrals &ints) {
  const auto Bk = ints.Bk.emptyQ() ? nullptr : &ints.Bk;
  const auto Sk = ints.Sk.emptyQ() ? nullptr : &ints.Sk;
  const auto dSk = ints.dSk.emptyQ() ? nullptr : &ints.dSk;
  const auto Sigma1_in_CI =
    ints.Sigma1_in_CI.empty() ? nullptr : &ints.Sigma1_in_CI;
  // No E0 needed: no dSigma/dE correction, and no BW Sigma_2
  if (ints.s1_corr.empty() && dSk == nullptr) {
    return construct_Hci(psi, ints.h1, ints.qk, Bk, Sk, nullptr, ints.hk,
                         nullptr, 0.0, Sigma1_in_CI);
  }
  // E0 is per (J, parity): if psi has been solved, it is just the lowest
  // energy of this J/Pi. Otherwise, fall back to the stored E0.
  auto s1_Jpi = ints.s1_corr;
  if (psi.num_solutions() > 0) {
    s1_Jpi.E0 = psi.energy(0);
  }
  return construct_Hci(psi, ints.h1, ints.qk, Bk, Sk, &s1_Jpi, ints.hk, dSk,
                       s1_Jpi.E0 - ints.E0_sigma2, Sigma1_in_CI);
}

} // namespace CI