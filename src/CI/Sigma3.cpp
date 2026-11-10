#include "Sigma3.hpp"
#include "Angular/include.hpp"
#include "CI_Integrals.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace CI {

//==============================================================================

Sigma3Orbitals sigma3_orbitals(const std::vector<DiracSpinor> &core,
                               const std::vector<DiracSpinor> &ci_basis) {
  Sigma3Orbitals out;

  std::size_t max_index = 0;
  for (const auto &a : core) {
    out.core.push_back(a.nk_index());
    max_index = std::max(max_index, std::size_t(a.nk_index()));
  }
  for (const auto &v : ci_basis) {
    max_index = std::max(max_index, std::size_t(v.nk_index()));
  }
  out.en.assign(max_index + 1, 0.0);
  out.en_bar.assign(max_index + 1, 0.0);

  for (const auto &a : core) {
    out.en[a.nk_index()] = a.en();
  }
  for (const auto &v : ci_basis) {
    out.en[v.nk_index()] = v.en();
    // Lowest CI orbital of this kappa
    auto e_bar = v.en();
    for (const auto &w : ci_basis) {
      if (w.kappa() == v.kappa() && w.en() < e_bar) {
        e_bar = w.en();
      }
    }
    out.en_bar[v.nk_index()] = e_bar;
  }
  out.en_min = ci_basis.empty() ? 0.0 : DiracSpinor::min_En(ci_basis);

  return out;
}

//==============================================================================

double Sigma3_inverse_denominator(MBPT::Denominators denominators,
                                  const Sigma3Orbitals &orbitals,
                                  DiracSpinor::Index a, DiracSpinor::Index v,
                                  DiracSpinor::Index w, DiracSpinor::Index u,
                                  DiracSpinor::Index x, DiracSpinor::Index y,
                                  DiracSpinor::Index z, double E0) {
  const auto &en = orbitals.en;
  const auto &en_bar = orbitals.en_bar;
  const auto ea = en[a];

  switch (denominators) {
  case MBPT::Denominators::RS:
    return 0.5 * (1.0 / (ea + en[u] - en[y] - en[z]) +
                  1.0 / (ea + en[x] - en[v] - en[w]));
  case MBPT::Denominators::Fermi:
    return 0.5 * (1.0 / (ea + en_bar[u] - en_bar[y] - en_bar[z]) +
                  1.0 / (ea + en_bar[x] - en_bar[v] - en_bar[w]));
  case MBPT::Denominators::Fermi0:
    return 1.0 / (ea - orbitals.en_min);
  case MBPT::Denominators::DFK: {
    const auto e_int = en[v] + en[w] + en[y] + en[z] - ea;
    return 0.5 * (1.0 / (en_bar[v] + en_bar[w] + en_bar[u] - e_int) +
                  1.0 / (en_bar[y] + en_bar[z] + en_bar[x] - e_int));
  }
  case MBPT::Denominators::BW:
    return 1.0 / (E0 - (en[v] + en[w] + en[y] + en[z] - ea));
  }
  assert(false);
  return 0.0;
}

//==============================================================================

double Sigma3_parents(const CSF3::Parent &bra, const CSF3::Parent &ket,
                      int twoJ, const Coulomb::QkTable &qk,
                      const Sigma3Orbitals &orbitals,
                      MBPT::Denominators denominators, double E0) {

  // Ket: pair (v, w) with J2, and u. Bra: x, and pair (y, z) with K
  const auto [v, w] = ket.pair;
  const auto u = ket.spectator;
  const auto twoJ2 = ket.twoJ2;
  const auto [y, z] = bra.pair;
  const auto x = bra.spectator;
  const auto twoK = bra.twoJ2;

  const auto tjx = Angular::nkindex_to_twoj(x);
  const auto tju = Angular::nkindex_to_twoj(u);

  double sum = 0.0;
  for (const auto a : orbitals.core) {
    const auto tja = Angular::nkindex_to_twoj(a);

    const auto sixj = Angular::sixj_2(tjx, tja, twoJ2, tju, twoJ, twoK);
    if (sixj == 0.0)
      continue;

    const auto g_yzau = CSF2_Coulomb(qk, y, z, a, u, twoK);
    if (g_yzau == 0.0)
      continue;
    const auto g_xavw = CSF2_Coulomb(qk, x, a, v, w, twoJ2);
    if (g_xavw == 0.0)
      continue;

    const auto inv_de = Sigma3_inverse_denominator(denominators, orbitals, a, v,
                                                   w, u, x, y, z, E0);

    sum +=
      Angular::neg1pow_2(twoK + tja + tju) * sixj * g_yzau * g_xavw * inv_de;
  }

  // 1/2: the pairs (x, a) and (a, u) are not antisymmetrised
  return 0.5 * std::sqrt((twoJ2 + 1.0) * (twoK + 1.0)) * sum;
}

//==============================================================================

double Sigma3_AB(const CSF3 &A, const CSF3 &B, int twoJ,
                 const Coulomb::QkTable &qk, const Sigma3Orbitals &orbitals,
                 MBPT::Denominators denominators, double E0) {
  double out = 0.0;
  for (const auto &p : A.parents()) {
    for (const auto &q : B.parents()) {
      out += p.coef * q.coef *
             Sigma3_parents(p, q, twoJ, qk, orbitals, denominators, E0);
    }
  }
  // 3!: orderings of the three electrons among the roles in the diagram
  return 6.0 * out;
}

} // namespace CI
