#include "Amplitudes/SecondOrder.hpp"
#include "Angular/include.hpp"
#include "CI/SecondOrder.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "ExternalField/TDHF.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/ostream.hpp"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace Amplitudes {

//==============================================================================
bool allowed_K(int K, int kt, int ks, int twoJb, int twoJa) {
  return K >= std::abs(kt - ks) && K <= kt + ks &&
         Angular::triangle(twoJb, 2 * K, twoJa) != 0;
}

//==============================================================================
int smallest_allowed_K(int kt, int ks, int twoJb, int twoJa) {
  for (int K = std::abs(kt - ks); K <= kt + ks; ++K) {
    if (allowed_K(K, kt, ks, twoJb, twoJa))
      return K;
  }
  return -1;
}

//==============================================================================
// ME from the table if present, otherwise calculated directly (with RPA; conj
// selects the dV component, see CorePolarisation::dV)
static double table_me(const Coulomb::meTable<double> &tab,
                       const DiracOperator::TensorOperator *h,
                       const ExternalField::CorePolarisation *dV,
                       const DiracSpinor &x, const DiracSpinor &y, bool conj) {
  const auto *ptr = tab.empty() ? nullptr : tab.get(x, y);
  if (ptr)
    return *ptr;
  return h->reducedME(x, y) + (dV ? dV->dV(x, y, conj) : 0.0);
}

//==============================================================================
double sos_valence(int K, const DiracSpinor &Fb, const DiracSpinor &Fa,
                   const DiracOperator::TensorOperator *t,
                   const DiracOperator::TensorOperator *s, double omega,
                   double omega_s, const std::vector<DiracSpinor> &spectrum,
                   const ExternalField::CorePolarisation *dVt,
                   const ExternalField::CorePolarisation *dVs,
                   const Coulomb::meTable<double> &t_me,
                   const Coulomb::meTable<double> &s_me) {

  if (Fb.parity() * Fa.parity() != t->parity() * s->parity())
    return 0.0;
  if (!allowed_K(K, t->rank(), s->rank(), Fb.twoj(), Fa.twoj()))
    return 0.0;

  const auto kt = t->rank();
  const auto ks = s->rank();

  // Each vertex carries the RPA component at the frequency of its operator:
  // the emission one for a negative frequency
  const auto conj_t = omega < 0.0;
  const auto conj_s = omega_s < 0.0;

  double A = 0.0;

  for (const auto &n : spectrum) {

    const auto [c1, c2] =
      CI::A_K_coefs(K, kt, ks, Fb.twoj(), n.twoj(), Fa.twoj());

    // 'ts' term: s takes a -> n, t takes n -> b
    if (c1 != 0.0 && !t->isZero(Fb, n) && !s->isZero(n, Fa)) {
      const auto denom = Fa.en() + omega_s - n.en();
      A += c1 * table_me(t_me, t, dVt, Fb, n, conj_t) *
           table_me(s_me, s, dVs, n, Fa, conj_s) / denom;
    }

    // 'st' term: t takes a -> n, s takes n -> b
    if (c2 != 0.0 && !s->isZero(Fb, n) && !t->isZero(n, Fa)) {
      const auto denom = Fa.en() + omega - n.en();
      A += c2 * table_me(s_me, s, dVs, Fb, n, conj_s) *
           table_me(t_me, t, dVt, n, Fa, conj_t) / denom;
    }
  }

  return A;
}

//==============================================================================
double sos_core(int K, int twoJ, const DiracOperator::TensorOperator *t,
                const DiracOperator::TensorOperator *s, double omega,
                double omega_s, const std::vector<DiracSpinor> &core,
                const std::vector<DiracSpinor> &excited,
                const ExternalField::CorePolarisation *dVt,
                const ExternalField::CorePolarisation *dVs) {
  // Identical quantity; single-particle j plays the role of J
  return CI::A_K_core(K, twoJ, t, s, omega, omega_s, core, excited, dVt, dVs);
}

//==============================================================================
// Projects each mixed state onto the span of the given states; e.g., the HF
// core, which selects the core-valence (Pauli blocking) part of the sum
static std::vector<DiracSpinor>
project_span(std::vector<DiracSpinor> dPsis,
             const std::vector<DiracSpinor> &states) {
  if (states.empty())
    return dPsis;
  for (auto &x : dPsis) {
    auto x_proj = 0.0 * x;
    for (const auto &c : states) {
      if (c.kappa() == x.kappa()) {
        x_proj += (c * x) * c;
      }
    }
    x = x_proj;
  }
  return dPsis;
}

//==============================================================================
std::pair<double, double>
ms_valence(int K, const DiracSpinor &Fb, const DiracSpinor &Fa,
           const DiracOperator::TensorOperator *t,
           const DiracOperator::TensorOperator *s, double omega, double omega_s,
           const ExternalField::TDHF *dVt, const ExternalField::TDHF *dVs,
           const MBPT::CorrelationPotential *Sigma,
           const std::vector<DiracSpinor> &project_onto,
           std::ostream &outstream) {

  assert(dVt != nullptr && dVs != nullptr &&
         "ms_valence requires TDHF objects for both operators");

  if (Fb.parity() * Fa.parity() != t->parity() * s->parity())
    return {0.0, 0.0};
  if (!allowed_K(K, t->rank(), s->rank(), Fb.twoj(), Fa.twoj()))
    return {0.0, 0.0};

  const auto kt = t->rank();
  const auto ks = s->rank();
  const auto twoJa = Fa.twoj();
  const auto twoJb = Fb.twoj();
  using ExternalField::dPsiType;

  // Sanity: energy conservation. A small violation (e.g. omega from
  // experiment) shifts the denominators of the two routes differently
  const auto de = std::abs(Fb.en() - Fa.en() - omega - omega_s);
  if (de > 1.0e-10) {
    fmt::print(outstream,
               "Warning: ms_valence: energy conservation violated by {:.2e}: "
               "omega + omega_s should be e_b - e_a\n",
               de);
  }

  // Each mixed state below is an X-type solution: for state v and operator h
  // at frequency w, dv = sum_n |n> <n||h + dV_h||v> / (e_v + w - e_n),
  // decomposed into kappa channels (one spinor per channel)

  // Every vertex carries the RPA component at the frequency of its operator
  // (the emission one for a negative frequency). The outer vertices take it
  // explicitly. The inner ones have it already: solve_dPsi takes dV^dag at a
  // negative frequency, and the symmetry sign on the mixed states of b
  // conjugates once more, so with energy conservation each inner vertex is
  // dV at its operator's frequency
  const auto conj_t = omega < 0.0;
  const auto conj_s = omega_s < 0.0;

  // Route 1: the sums formed with the mixed states of s
  double A_s = 0.0;
  {
    // 'ts' term: da = mixed state of a due to s at omega_s:
    // denominators e_a + omega_s - e_n
    const auto da = project_span(
      dVs->solve_dPsis(Fa, omega_s, dPsiType::X, Sigma), project_onto);
    for (const auto &da_x : da) {
      const auto [c1, c2] = CI::A_K_coefs(K, kt, ks, twoJb, da_x.twoj(), twoJa);
      if (c1 != 0.0 && !t->isZero(Fb, da_x)) {
        A_s += c1 * (t->reducedME(Fb, da_x) + dVt->dV(Fb, da_x, conj_t));
      }
    }
    // 'st' term: db = mixed state of b due to s at e_a + omega - e_b:
    // denominators e_a + omega - e_n. The solve gives <n||s+dVs||b>; the
    // symmetry sign converts to the required <b||s+dVs||n>
    const auto db = project_span(
      dVs->solve_dPsis(Fb, Fa.en() + omega - Fb.en(), dPsiType::X, Sigma),
      project_onto);
    for (const auto &db_x : db) {
      const auto [c1, c2] = CI::A_K_coefs(K, kt, ks, twoJb, db_x.twoj(), twoJa);
      if (c2 != 0.0 && !t->isZero(db_x, Fa)) {
        A_s += c2 * CI::symm_sign(s, twoJb, db_x.twoj()) *
               (t->reducedME(db_x, Fa) + dVt->dV(db_x, Fa, conj_t));
      }
    }
  }

  // Route 2: the same sums, formed with the mixed states of t
  double A_t = 0.0;
  {
    // 'st' term: Da = mixed state of a due to t at omega:
    // denominators e_a + omega - e_n
    const auto Da = project_span(
      dVt->solve_dPsis(Fa, omega, dPsiType::X, Sigma), project_onto);
    for (const auto &Da_x : Da) {
      const auto [c1, c2] = CI::A_K_coefs(K, kt, ks, twoJb, Da_x.twoj(), twoJa);
      if (c2 != 0.0 && !s->isZero(Fb, Da_x)) {
        A_t += c2 * (s->reducedME(Fb, Da_x) + dVs->dV(Fb, Da_x, conj_s));
      }
    }
    // 'ts' term: Db = mixed state of b due to t at e_a + omega_s - e_b:
    // denominators e_a + omega_s - e_n
    const auto Db = project_span(
      dVt->solve_dPsis(Fb, Fa.en() + omega_s - Fb.en(), dPsiType::X, Sigma),
      project_onto);
    for (const auto &Db_x : Db) {
      const auto [c1, c2] = CI::A_K_coefs(K, kt, ks, twoJb, Db_x.twoj(), twoJa);
      if (c1 != 0.0 && !s->isZero(Db_x, Fa)) {
        A_t += c1 * CI::symm_sign(t, twoJb, Db_x.twoj()) *
               (s->reducedME(Db_x, Fa) + dVs->dV(Db_x, Fa, conj_s));
      }
    }
  }

  return {A_s, A_t};
}

//==============================================================================
double ms_core(int K, int twoJ, const DiracOperator::TensorOperator *t,
               const DiracOperator::TensorOperator *s, double omega,
               double omega_s, const std::vector<DiracSpinor> &core,
               const ExternalField::TDHF *dVt, const ExternalField::TDHF *dVs,
               const MBPT::CorrelationPotential *Sigma) {

  assert(dVt != nullptr && dVs != nullptr &&
         "ms_core requires TDHF objects for both operators");

  // The core is a closed shell: only a scalar, even-parity, amplitude survives
  if (K != 0 || t->parity() * s->parity() != 1)
    return 0.0;
  if (core.empty())
    return 0.0;

  const auto kt = t->rank();
  const auto ks = s->rank();
  using ExternalField::dPsiType;

  double A_core = 0.0;
  for (const auto &Fc : core) {
    double Uc = 0.0;

    // 'ts': s excites the core electron (mixed state carries dVs); t, bare,
    // returns it. Denominators e_c + omega_s - e_m
    const auto Xs = dVs->solve_dPsis(Fc, omega_s, dPsiType::X, Sigma);
    for (const auto &x : Xs) {
      const auto [c1, c2] =
        CI::A_K_coefs(0, kt, ks, Fc.twoj(), x.twoj(), Fc.twoj());
      if (c1 != 0.0 && !t->isZero(Fc, x)) {
        Uc += c1 * t->reducedME(Fc, x);
      }
    }

    // 'st': t excites the core electron (mixed state carries dVt); s returns
    // it. Denominators e_c + omega - e_m
    const auto Xt = dVt->solve_dPsis(Fc, omega, dPsiType::X, Sigma);
    for (const auto &x : Xt) {
      const auto [c1, c2] =
        CI::A_K_coefs(0, kt, ks, Fc.twoj(), x.twoj(), Fc.twoj());
      if (c2 != 0.0 && !s->isZero(Fc, x)) {
        Uc += c2 * s->reducedME(Fc, x);
      }
    }

    A_core += std::sqrt(double(Fc.twojp1())) * Uc;
  }

  return std::sqrt(double(twoJ + 1)) * A_core;
}

} // namespace Amplitudes
