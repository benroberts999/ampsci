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
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
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
double CSF2_Sigma2(const Coulomb::LkTable &Sk, DiracSpinor::Index v,
                   DiracSpinor::Index w, DiracSpinor::Index x,
                   DiracSpinor::Index y, int twoJ, const Coulomb::QkTable *qk,
                   const std::vector<double> &hk, const Coulomb::LkTable *dSk,
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
                 const Coulomb::LkTable &Sk, const Coulomb::QkTable *qk,
                 const std::vector<double> &hk, const Coulomb::LkTable *dSk,
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
double Sigma1Correction::delta_h1(DiracSpinor::Index a, DiracSpinor::Index b,
                                  DiracSpinor::Index spectator) const {
  const auto Sig = S1.getv(a, b);
  // E0 = 0.0 means it was never set: no correction (see iterate_E0)
  if (Sig == 0.0 || E0 == 0.0) {
    return 0.0;
  }
  const auto es = e_sigma.find(Angular::nkindex_to_kappa(a));
  const auto e_spec = en.find(spectator);
  if (es == e_sigma.end() || e_spec == en.end()) {
    return 0.0;
  }
  const auto dE = E0 - es->second - e_spec->second;
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
LinAlg::Matrix<double>
iterate_E0(PsiJPi *psi, const Sigma1Correction &s1c,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
           const Coulomb::WkTable *Bk, const Coulomb::LkTable *Sk,
           const std::vector<double> &hk, const Coulomb::LkTable *dSk,
           double E0_sigma2, std::ostream &outstream) {
  assert(psi != nullptr);

  constexpr auto eps_E0 = 1.0e-8;
  constexpr auto max_iterations = 16;
  const auto N_CSFs = psi->CSFs().size();

  // Sigma_2 was tabulated at the reference E0_sigma2, and is shifted to the
  // current E0 using dS^k/dE0 (see corrected_Sk). Both corrections then move
  // together as E0 is iterated.
  const auto dE0_s2 = [dSk, E0_sigma2](double E0) {
    return dSk == nullptr ? 0.0 : E0 - E0_sigma2;
  };

  // E0 belongs to this (J, parity), but the correction tables are shared
  // between them: work on a local copy. Only E0 changes between passes.
  auto s1_block = s1c;

  // E0 = 0.0 means "not set": start from the lowest zeroth-order configuration
  // (two electrons in the lowest orbital)
  if (s1_block.E0 == 0.0) {
    auto e_min = 0.0;
    for (const auto &orbital : s1_block.en) {
      if (orbital.second < e_min) {
        e_min = orbital.second;
      }
    }
    s1_block.E0 = 2.0 * e_min;
  }

  fmt::print(outstream,
             "Iterating E0 for dSigma/dE correction (initial: {:.8f} au):\n",
             s1_block.E0);
  outstream << std::flush;

  // First pass: build Hci at the current E0, and diagonalise for the lowest
  // level. E0 enters only through the (small) dSigma/dE correction, so the
  // state itself hardly changes from pass to pass. Later passes therefore
  // rebuild Hci with the updated E0, and take the new E0 as the expectation
  // value of the new Hci in the (unchanged) state: first-order perturbation
  // theory in the change to Hci. Building Hci is N^2; diagonalising is N^3,
  // so this is much cheaper than re-solving each pass.
  auto Hci = construct_Hci(*psi, h1, qk, Bk, Sk, &s1_block, hk, dSk,
                           dE0_s2(s1_block.E0));
  psi->solve(Hci, 1);
  const auto c = psi->coefs(0);

  auto delta = psi->energy(0) - s1_block.E0;
  s1_block.E0 = psi->energy(0);
  fmt::print(outstream, "  it {:2}: E0 = {:.8f} au (eps = {:.1e})\n", 1,
             s1_block.E0, delta / s1_block.E0);
  outstream << std::flush;

  for (int it = 2; it <= max_iterations && std::abs(delta) >= eps_E0; ++it) {
    Hci = construct_Hci(*psi, h1, qk, Bk, Sk, &s1_block, hk, dSk,
                        dE0_s2(s1_block.E0));

    // <c|Hci|c>; c is normalised
    double E0_new = 0.0;
#pragma omp parallel for reduction(+ : E0_new)
    for (std::size_t i = 0; i < N_CSFs; ++i) {
      double row = 0.0;
      for (std::size_t j = 0; j < N_CSFs; ++j) {
        row += Hci(i, j) * c[j];
      }
      E0_new += c[i] * row;
    }

    delta = E0_new - s1_block.E0;
    s1_block.E0 = E0_new;
    fmt::print(outstream, "  it {:2}: E0 = {:.8f} au (eps = {:.1e})\n", it,
               s1_block.E0, delta / s1_block.E0);
    outstream << std::flush;
    if (it == max_iterations && std::abs(delta) >= eps_E0) {
      outstream << "Warning: E0 iteration did not converge\n";
    }
  }

  return Hci;
}

//==============================================================================
// Determines CI Hamiltonian matrix element for two 2-particle CSFs, a and b
double Hab(const CI::CSF2 &X, const CI::CSF2 &V, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
           const Sigma1Correction *s1c) {

  // Calculates matrix element of the CI Hamiltonian between two CSFs

  const auto [v, w] = V.states;
  const auto [x, y] = X.states;

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

  return h1_VX + CSF2_Coulomb(qk, v, w, x, y, twoJ);
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
                 const std::vector<CI::CSF2> &CSFAs, int twoJA,
                 const LinAlg::View<const double> &cB,
                 const std::vector<CI::CSF2> &CSFBs, int twoJB,
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
  const auto NA = CSFAs.size();
  const auto NB = CSFBs.size();

  double rme = 0.0;
#pragma omp parallel for collapse(2) reduction(+ : rme)
  for (std::size_t i = 0ul; i < NA; ++i) {
    for (std::size_t j = 0ul; j < NB; ++j) {
      const auto &csf_i = CSFAs.at(i);
      const auto ci = cA[i];
      const auto &csf_j = CSFBs.at(j);
      const auto cj = cB[j];

      rme += ci * cj * RME_CSF2(csf_i, twoJA, csf_j, twoJB, h, K_rank);
    }
  }
  return rme;
}

//==============================================================================
// Calculate reduce ME between two 2-particle CSFs
double RME_CSF2(const CI::CSF2 &X, int twoJX, const CI::CSF2 &V, int twoJV,
                const Coulomb::meTable<double> &h, int K_rank) {

  const auto [v, w] = V.states;
  const auto [x, y] = X.states;
  const auto etaV = v == w ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto etaX = x == y ? 1.0 / std::sqrt(2.0) : 1.0;

  // const auto num_diff = CI::CSF2::num_different(F, I);
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
std::pair<std::string, double>
leading_config(const LinAlg::View<const double> &coefs,
               const std::vector<CSF2> &csfs) {
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
  // Nearest (S, L) consistent with J; for two electrons, S = 0 or 1.
  // L limit: generous; l1 + l2 for any realistic CI basis
  constexpr int max_L = 16;

  int best_L = -1;
  int best_S = 0;
  double best_del = 1.0e30;
  for (int S = 0; S <= 1; ++S) {
    for (int L = 0; L <= max_L; ++L) {
      if (twoJ > 2 * (L + S) || twoJ < 2 * std::abs(L - S))
        continue;
      const auto del =
        std::abs(L * (L + 1.0) - L2) + std::abs(S * (S + 1.0) - S2);
      if (del < best_del) {
        best_del = del;
        best_L = L;
        best_S = S;
      }
    }
  }
  return {best_S, best_L};
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
LinAlg::Matrix<double>
construct_Hci(const PsiJPi &psi, const Coulomb::meTable<double> &h1,
              const Coulomb::QkTable &qk, const Coulomb::WkTable *Bk,
              const Coulomb::LkTable *Sk, const Sigma1Correction *s1c,
              const std::vector<double> &hk, const Coulomb::LkTable *dSk,
              double dE0) {

  const auto N_CSFs = psi.CSFs().size();
  const auto twoJ = psi.twoJ();
  LinAlg::Matrix Hci(N_CSFs, N_CSFs);

#pragma omp parallel for
  for (std::size_t iA = 0; iA < N_CSFs; ++iA) {
    // go to iB <= iA only: symmetric matrix
    for (std::size_t iB = 0; iB <= iA; ++iB) {

      const auto &A = psi.CSF(iA);
      const auto &B = psi.CSF(iB);

      // Regular CI matrix (h1 may include Sigma_1):
      const auto E_AB = Hab(A, B, twoJ, h1, qk, s1c);
      // Sigma_2 correction:
      const auto dEs_AB =
        Sk ? CI::Sigma2_AB(A, B, twoJ, *Sk, &qk, hk, dSk, dE0) : 0.0;
      // Breit correction:
      const auto dEb_AB = Bk ? CI::Breit_AB(A, B, twoJ, *Bk) : 0.0;

      Hci(iA, iB) = E_AB + dEs_AB + dEb_AB;
      // fill other half of symmetric matrix:
      if (iB != iA) {
        Hci(iB, iA) = E_AB + dEs_AB + dEb_AB;
      }
    }
  }

  return Hci;
}

//==============================================================================
LinAlg::Matrix<double> construct_Hci(const PsiJPi &psi, const Integrals &ints) {
  const auto Bk = ints.Bk.emptyQ() ? nullptr : &ints.Bk;
  const auto Sk = ints.Sk.emptyQ() ? nullptr : &ints.Sk;
  const auto dSk = ints.dSk.emptyQ() ? nullptr : &ints.dSk;
  if (ints.s1_corr.empty()) {
    return construct_Hci(psi, ints.h1, ints.qk, Bk, Sk, nullptr, ints.hk);
  }
  // E0 is per (J, parity): if psi has been solved, it is just the lowest
  // energy of this block. Otherwise, fall back to the stored E0.
  auto s1_block = ints.s1_corr;
  if (psi.num_solutions() > 0) {
    s1_block.E0 = psi.energy(0);
  }
  return construct_Hci(psi, ints.h1, ints.qk, Bk, Sk, &s1_block, ints.hk, dSk,
                       s1_block.E0 - ints.E0_sigma2);
}

} // namespace CI