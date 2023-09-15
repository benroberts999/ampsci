#include "CI_Integrals.hpp"
#include "CSF.hpp"
#include "Coulomb/Coulomb.hpp"
#include "MBPT/Sigma2.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>

namespace CI {

//==============================================================================
// Calculates the anti-symmetrised Coulomb integral for 2-particle states:
// C1*C2*(g_abcd-g_abdc), where Cs are C.G. coefficients
double CSF2_Coulomb(const Coulomb::QkTable &qk, const DiracSpinor &v,
                    const DiracSpinor &w, const DiracSpinor &x,
                    const DiracSpinor &y, int twoJ) {

  // If c==d, or a==b : can make short-cut due to symmetry
  // More efficient to use two Q's than W:

  double out = 0.0;

  // Direct part:
  const auto [k0, k1] = Coulomb::k_minmax_Q(v, w, x, y);
  for (int k = k0; k <= k1; k += 2) {
    const auto sjs =
        Angular::sixj_2(v.twoj(), w.twoj(), twoJ, y.twoj(), x.twoj(), 2 * k);
    if (sjs == 0.0)
      continue;
    const auto qk_abcd = qk.Q(k, v, w, x, y);
    const auto s = Angular::neg1pow_2(v.twoj() + x.twoj() + 2 * k + twoJ);
    out += s * sjs * qk_abcd;
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
  const auto [l0, l1] = Coulomb::k_minmax_Q(v, w, y, x);
  for (int k = l0; k <= l1; k += 2) {
    const auto sjs =
        Angular::sixj_2(v.twoj(), w.twoj(), twoJ, x.twoj(), y.twoj(), 2 * k);
    if (sjs == 0.0)
      continue;
    const auto qk_abdc = qk.Q(k, v, w, y, x);
    const auto s = Angular::neg1pow_2(v.twoj() + x.twoj() + 2 * k);
    out += s * sjs * qk_abdc;
  }

  return out;
}

//==============================================================================
double CSF2_Sigma2(const Coulomb::LkTable &Sk, const DiracSpinor &v,
                   const DiracSpinor &w, const DiracSpinor &x,
                   const DiracSpinor &y, int twoJ) {

  // If c==d, or a==b : can make short-cut due to symmetry
  // More efficient to use two Q's than W:

  double out = 0.0;

  // Direct part:
  const auto [k0, k1] = MBPT::k_minmax_S(v, w, x, y);
  for (int k = k0; k <= k1; ++k) {
    const auto sjs =
        Angular::sixj_2(v.twoj(), w.twoj(), twoJ, y.twoj(), x.twoj(), 2 * k);
    if (sjs == 0.0)
      continue;
    // const auto Sk_abcd = MBPT::Sk_vwxy(k, a, b, c, d, qk, core, excited, SixJ);
    const auto Sk_abcd = Sk.Q(k, v, w, x, y);
    const auto s = Angular::neg1pow_2(v.twoj() + x.twoj() + 2 * k + twoJ);
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
  const auto [l0, l1] = MBPT::k_minmax_S(v, w, y, x);
  for (int k = l0; k <= l1; ++k) {
    const auto sjs =
        Angular::sixj_2(v.twoj(), w.twoj(), twoJ, x.twoj(), y.twoj(), 2 * k);
    if (sjs == 0.0)
      continue;
    // const auto Sk_abdc = MBPT::Sk_vwxy(k, a, b, d, c, qk, core, excited, SixJ);
    const auto Sk_abdc = Sk.Q(k, v, w, y, x);
    const auto s = Angular::neg1pow_2(v.twoj() + x.twoj() + 2 * k);
    out += s * sjs * Sk_abdc;
  }

  return out;
}

//==============================================================================
double Sigma2_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                 const Coulomb::LkTable &Sk) {
  const auto [v, w] = A.states;
  const auto [x, y] = B.states;
  return CSF2_Sigma2(Sk, *v, *w, *x, *y, twoJ);
}

//==============================================================================
// Determines CI Hamiltonian matrix element for two 2-particle CSFs, a and b
double Hab(const CI::CSF2 &V, const CI::CSF2 &X, int twoJ,
           const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk) {

  // Calculates matrix element of the CI Hamiltonian between two CSFs

  const auto [v, w] = V.states;
  const auto [x, y] = X.states;

  const auto etaV = v == w ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto etaX = x == y ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto eta2 = etaV * etaX;

  double h1_VX = 0.0;
  if (y == w) {
    h1_VX += eta2 * h1.getv(*v, *x);
  }
  if (y == v) {
    h1_VX += eta2 * h1.getv(*w, *x);
  }
  if (x == w) {
    h1_VX += eta2 * h1.getv(*v, *y);
  }
  if (x == v) {
    h1_VX += eta2 * h1.getv(*w, *y);
  }

  // understand why these two are different => the answer!
  // BUT, they are the same when no Sigma included..

  // double h1_VX = 0.0;
  // auto nd = CI::CSF2::num_different(V, X);
  // if (nd == 0) {
  //   h1_VX = h1.getv(*v, *v) + h1.getv(*w, *w);
  // } else if (nd == 1) {
  //   auto [n, a] = CI::CSF2::diff_1_na(V, X);
  //   assert(n != a);
  //   h1_VX = eta2 * h1.getv(*n, *a);
  // }

  return h1_VX + CSF2_Coulomb(qk, *v, *w, *x, *y, twoJ);
}

//==============================================================================
Coulomb::meTable<double>
calculate_h1_table(const std::vector<DiracSpinor> &ci_basis,
                   const std::vector<DiracSpinor> &s1_basis_core,
                   const std::vector<DiracSpinor> &s1_basis_excited,
                   const Coulomb::QkTable &qk, bool include_Sigma1) {
  // Create lookup table for one-particle matrix elements, h1
  Coulomb::meTable<double> h1;

  // Calculate + store all 1-body integrals
  for (const auto &v : ci_basis) {

    // Find lowest valence state of this l; for Sigma energy
    const auto vp =
        std::find_if(ci_basis.begin(), ci_basis.end(),
                     [l = v.l()](const auto &n) { return n.l() == l; });
    auto ev = vp->en();

    for (const auto &w : ci_basis) {
      if (w > v)
        continue;
      if (w.kappa() != v.kappa())
        continue;
      const auto h0_vw = v == w ? v.en() : 0.0;

      // Can use Sigma matrix instead: all-orders?
      const auto Sigma_vw = include_Sigma1 ?
                                MBPT::Sigma_vw(v, w, qk, s1_basis_core,
                                               s1_basis_excited, 99, ev) :
                                0.0;

      h1.add(v, w, h0_vw + Sigma_vw);
      // Add symmetric partner:
      if (v != w)
        h1.add(w, v, h0_vw + Sigma_vw);
    }
  }
  return h1;
}

//==============================================================================
Coulomb::LkTable calculate_Sk(const std::string &filename,
                              const std::vector<DiracSpinor> &cis2_basis,
                              const std::vector<DiracSpinor> &s2_basis_core,
                              const std::vector<DiracSpinor> &s2_basis_excited,
                              const Coulomb::QkTable &qk, int max_k,
                              bool exclude_wrong_parity_box) {

  Coulomb::LkTable Sk;

  const auto max_twoj = std::max({DiracSpinor::max_tj(s2_basis_excited),
                                  DiracSpinor::max_tj(s2_basis_core),
                                  DiracSpinor::max_tj(cis2_basis)});
  Angular::SixJTable sjt(max_twoj);

  const auto Sk_function = [&](int k, const DiracSpinor &v,
                               const DiracSpinor &w, const DiracSpinor &x,
                               const DiracSpinor &y) {
    return MBPT::Sk_vwxy(k, v, w, x, y, qk, s2_basis_core, s2_basis_excited,
                         sjt);
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
    Sk.fill(cis2_basis, Sk_function, Sk_selection_rule, max_k);

    const auto total = Sk.count();
    assert(total >= existing);
    const auto new_integrals = total - existing;
    std::cout << "Calculated " << new_integrals << " new MBPT integrals\n";
    if (new_integrals > 0) {
      Sk.write(filename);
    }
  }
  return Sk;
}

//==============================================================================
// Takes a subset of input basis according to subset_string.
// Only states *not* included in frozen_core_string are included.
std::vector<DiracSpinor> basis_subset(const std::vector<DiracSpinor> &basis,
                                      const std::string &subset_string,
                                      const std::string &frozen_core_string) {

  // Form 'subset' from {a} in 'basis', if:
  //    a _is_ in subset_string AND
  //    a _is not_ in basis string

  std::vector<DiracSpinor> subset;
  const auto nmaxk_list = AtomData::n_kappa_list(subset_string);
  const auto core_list = AtomData::core_parser(frozen_core_string);

  for (const auto &a : basis) {

    // Check if a is present in 'subset_string'
    const auto nk =
        std::find_if(nmaxk_list.cbegin(), nmaxk_list.cend(),
                     [&a](const auto &tnk) { return a.kappa() == tnk.second; });
    if (nk == nmaxk_list.cend())
      continue;
    // nk is now max n, for given kappa {max_n, kappa}
    if (a.n() > nk->first)
      continue;

    // assume only filled shells in frozen core
    const auto core = std::find_if(
        core_list.cbegin(), core_list.cend(), [&a](const auto &tcore) {
          return a.n() == tcore.n && a.l() == tcore.l;
        });

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
                 const DiracOperator::TensorOperator *h) {

  // selection rules: Not required (operator encodes it's own),
  // but it's faster to check in advance
  const auto piA = CSFAs.front().parity();
  const auto piB = CSFBs.front().parity();
  if (piA * piB != h->parity())
    return 0.0;
  if (std::abs(twoJA - twoJB) > 2 * h->rank())
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

      rme += ci * cj * RME_CSF2(csf_i, twoJA, csf_j, twoJB, h);
    }
  }
  return rme;
}

//==============================================================================
// Overload for diagonal matrix elements
double ReducedME(const LinAlg::View<const double> &cA,
                 const std::vector<CI::CSF2> &CSFs, int twoJ,
                 const DiracOperator::TensorOperator *h) {
  return ReducedME(cA, CSFs, twoJ, cA, CSFs, twoJ, h);
}

//==============================================================================
// Calculate reduce ME between two 2-particle CSFs - XXX not quite right??
double RME_CSF2(const CI::CSF2 &V, int twoJV, const CI::CSF2 &X, int twoJX,
                const DiracOperator::TensorOperator *h) {

  const auto [v, w] = V.states;
  const auto [x, y] = X.states;
  const auto etaV = v == w ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto etaX = x == y ? 1.0 / std::sqrt(2.0) : 1.0;

  // const auto num_diff = CI::CSF2::num_different(V, X);
  const auto twok = 2 * h->rank();

  const auto f = etaV * etaX * std::sqrt(double(twoJV + 1) * (twoJX + 1)) *
                 Angular::neg1pow_2(twok);

  double sum = 0.0;
  if (y == w) {
    const auto sj =
        Angular::sixj_2(twoJV, twoJX, twok, v->twoj(), x->twoj(), w->twoj());
    const auto t = h->reducedME(*x, *v);
    const auto s = Angular::neg1pow_2(w->twoj() + x->twoj() + twoJV);
    sum += f * sj * t * s;
  }
  if (y == v) {
    const auto sj =
        Angular::sixj_2(twoJV, twoJX, twok, w->twoj(), x->twoj(), v->twoj());
    const auto t = h->reducedME(*x, *w);
    const auto s = Angular::neg1pow_2(w->twoj() + x->twoj());
    sum += f * sj * t * s;
  }
  if (x == w) {
    const auto sj =
        Angular::sixj_2(twoJV, twoJX, twok, v->twoj(), y->twoj(), w->twoj());
    const auto t = h->reducedME(*y, *v);
    const auto s = Angular::neg1pow_2(twoJV + twoJX + 1);
    sum += f * sj * t * s;
  }
  if (x == v) {
    const auto sj =
        Angular::sixj_2(twoJV, twoJX, twok, w->twoj(), y->twoj(), v->twoj());
    const auto t = h->reducedME(*y, *w); //yv?
    const auto s = Angular::neg1pow_2(w->twoj() + v->twoj() + twoJV);
    sum += f * sj * t * s;
  }
  return sum;
}

//==============================================================================
std::pair<int, int> Term_S_L(int l1, int l2, int twoJ, double gJ_target) {
  // Determine Term Symbol, from g-factor
  const auto min_L = std::abs(l1 - l2);
  const auto max_L = std::abs(l1 + l2);
  const auto min_S = 0;
  const auto max_S = 1;
  int L = min_L;
  int S = min_L;
  double best_del = 2.0;
  if (twoJ != 0) {
    for (int tL = min_L; tL <= max_L; ++tL) {
      for (int tS = min_S; tS <= max_S; ++tS) {
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
  // const auto tSp1 = 2.0 * S + 1.0;
  return {S, L};
}

} // namespace CI