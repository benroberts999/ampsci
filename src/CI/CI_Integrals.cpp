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
double CSF2_Coulomb(const Coulomb::QkTable &qk, const DiracSpinor &a,
                    const DiracSpinor &b, const DiracSpinor &c,
                    const DiracSpinor &d, int twoJ) {

  // If c==d, or a==b : can make short-cut due to symmetry
  // More efficient to use two Q's than W:

  double out = 0.0;

  // Direct part:
  const auto [k0, k1] = Coulomb::k_minmax_Q(a, b, c, d);
  for (int k = k0; k <= k1; k += 2) {
    const auto sjs =
        Angular::sixj_2(a.twoj(), b.twoj(), twoJ, d.twoj(), c.twoj(), 2 * k);
    if (sjs == 0.0)
      continue;
    const auto qk_abcd = qk.Q(k, a, b, c, d);
    const auto s = Angular::neg1pow_2(a.twoj() + c.twoj() + 2 * k + twoJ);
    out += s * sjs * qk_abcd;
  }

  // Take advantage of symmetries: faster (+ numerically stable)
  // c == d => J is even (identical states), eta2=1/sqrt(2)
  // eta_ab = 1/sqrt(2) if a==b
  // Therefore: e.g., if c==d
  // => eta_ab * eta_cd * (out + (-1)^J*out) = eta_ab * sqrt(2) * out
  if (a == b && c == d) {
    return out;
  } else if (c == d || a == b) {
    // by {ab},{cd} symmetry: same works for case a==b
    return std::sqrt(2.0) * out;
  }

  // Exchange part:
  const auto [l0, l1] = Coulomb::k_minmax_Q(a, b, d, c);
  for (int k = l0; k <= l1; k += 2) {
    const auto sjs =
        Angular::sixj_2(a.twoj(), b.twoj(), twoJ, c.twoj(), d.twoj(), 2 * k);
    if (sjs == 0.0)
      continue;
    const auto qk_abdc = qk.Q(k, a, b, d, c);
    const auto s = Angular::neg1pow_2(a.twoj() + c.twoj() + 2 * k);
    out += s * sjs * qk_abdc;
  }

  return out;
}

//==============================================================================
double CSF2_Sigma2(const Coulomb::QkTable &qk, const DiracSpinor &a,
                   const DiracSpinor &b, const DiracSpinor &c,
                   const DiracSpinor &d, int twoJ,
                   const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const Angular::SixJTable &SixJ) {

  // If c==d, or a==b : can make short-cut due to symmetry
  // More efficient to use two Q's than W:

  double out = 0.0;

  // Direct part:
  const auto [k0, k1] = Coulomb::k_minmax_Q(a, b, c, d);
  // for (int k = k0 ; k <= k1 ; += 2) { // No! diff rules!
  for (int k = 0; k <= k1 + 1; ++k) {
    const auto sjs =
        Angular::sixj_2(a.twoj(), b.twoj(), twoJ, d.twoj(), c.twoj(), 2 * k);
    if (sjs == 0.0)
      continue;
    const auto Sk_abcd = MBPT::Sk_vwxy(k, a, b, c, d, qk, core, excited, SixJ);
    const auto s = Angular::neg1pow_2(a.twoj() + c.twoj() + 2 * k + twoJ);
    out += s * sjs * Sk_abcd;
  }

  // Take advantage of symmetries: faster (+ numerically stable)
  // c == d => J is even (identical states), eta2=1/sqrt(2)
  // eta_ab = 1/sqrt(2) if a==b
  // Therefore: e.g., if c==d
  // => eta_ab * eta_cd * (out + (-1)^J*out) = eta_ab * sqrt(2) * out
  if (a == b && c == d) {
    return out;
  } else if (c == d || a == b) {
    // by {ab},{cd} symmetry: same works for case a==b
    return std::sqrt(2.0) * out;
  }

  // Exchange part:
  const auto [l0, l1] = Coulomb::k_minmax_Q(a, b, d, c);
  // for (int k = l0; k <= l1; k += 2) {
  for (int k = 0; k <= l1 + 1; ++k) {
    const auto sjs =
        Angular::sixj_2(a.twoj(), b.twoj(), twoJ, c.twoj(), d.twoj(), 2 * k);
    if (sjs == 0.0)
      continue;
    const auto Sk_abdc = MBPT::Sk_vwxy(k, a, b, d, c, qk, core, excited, SixJ);
    const auto s = Angular::neg1pow_2(a.twoj() + c.twoj() + 2 * k);
    out += s * sjs * Sk_abdc;
  }

  return out;
}

//==============================================================================
double Sigma2_AB(const CI::CSF2 &A, const CI::CSF2 &B, int twoJ,
                 const Coulomb::QkTable &qk,
                 const std::vector<DiracSpinor> &core,
                 const std::vector<DiracSpinor> &excited,
                 const Angular::SixJTable &SixJ) {
  const auto [v, w] = A.states;
  const auto [x, y] = B.states;
  return CSF2_Sigma2(qk, *v, *w, *x, *y, twoJ, core, excited, SixJ);
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
// Calculate reduced matrix elements between two CI states.
// cA is CI expansion coefficients (row if CI eigenvector matrix)
double ReducedME(const double *cA, const std::vector<CI::CSF2> &CSFAs,
                 int twoJA, const double *cB,
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
double ReducedME(const double *cA, const std::vector<CI::CSF2> &CSFs, int twoJ,
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

} // namespace CI