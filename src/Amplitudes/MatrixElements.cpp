#include "Amplitudes/MatrixElements.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "MBPT/StructureRad.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/ostream.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

namespace Amplitudes {

//==============================================================================
void solve_at_frequency(DiracOperator::TensorOperator *h,
                        ExternalField::CorePolarisation *dV,
                        MBPT::StructureRad *sr, double omega,
                        const MEoptions &options) {

  const bool op_each = options.operator_omega == Frequency::transition;
  const bool rpa_each = options.rpa_omega == Frequency::transition;

  // The RPA is solved with the operator at +|omega|
  if (op_each) {
    h->updateFrequency(std::abs(omega));
  }

  // Re-solve the RPA at this frequency; from scratch if the previous
  // solution is poorly converged, or if first-order RPA (which must not
  // iterate from a previous solution)
  if (rpa_each && dV) {
    if (dV->last_eps() > 1.0e-5 || options.rpa_iterations == 1 ||
        std::isnan(dV->last_eps())) {
      dV->clear();
    }
    dV->solve_core(omega, options.rpa_iterations, options.print);
  }

  // The element itself is taken with the operator at the signed frequency
  if (op_each) {
    h->updateFrequency(omega);
  }

  // The SR tables hold the operator and the RPA on the internal lines, so
  // they follow whatever was just changed
  if (sr && ((op_each && h->freqDependantQ()) || rpa_each)) {
    sr->solve_core(h, dV, omega);
  }
}

//==============================================================================
MEdata matrix_element(const DiracSpinor &a, const DiracSpinor &b,
                      const DiracOperator::TensorOperator *h,
                      const ExternalField::CorePolarisation *dV,
                      const MBPT::StructureRad *sr,
                      DiracOperator::MatrixElementType type, double omega,
                      double omega_sr) {

  MEdata me;
  me.a = a.shortSymbol();
  me.b = b.shortSymbol();
  me.omega = omega;
  me.factor = h->matel_factor(type, a, b);
  me.t0 = h->reducedME(a, b);
  // The RPA component at the transition frequency: emission if negative
  me.dv = dV ? dV->dV(a, b, omega < 0.0) : 0.0;
  me.has_rpa = (dV != nullptr);

  if (sr) {
    me.sr = sr->SR(a, b, omega_sr);
    me.norm_factor = sr->norm_factor(a, b);
    me.norm = me.norm_factor * (me.t0 + me.dv);
  }

  return me;
}

//==============================================================================
double sr_frequency(const ExternalField::CorePolarisation *dV, double omega,
                    const MEoptions &options) {
  const bool rpa_each = options.rpa_omega == Frequency::transition;
  return rpa_each ? omega : (dV ? dV->last_omega() : 0.0);
}

//==============================================================================
double frequency_derivative_term(const DiracSpinor &a, const DiracSpinor &b,
                                 DiracOperator::TensorOperator *h,
                                 const MBPT::StructureRad *sr, double omega) {
  if (!h->freqDependantQ() || a == b)
    return 0.0;

  // dT = (dt/dw) * dw^(2), with the derivative taken numerically
  const auto dw_2 = sr->Sigma_vw(a, a) - sr->Sigma_vw(b, b);
  const auto del = 0.05 * std::abs(omega);
  h->updateFrequency(omega + del);
  const auto t_plus = h->reducedME(a, b);
  h->updateFrequency(omega - del);
  const auto t_minus = h->reducedME(a, b);
  h->updateFrequency(omega);
  return (t_plus - t_minus) / (2 * del) * dw_2;
}

//==============================================================================
std::vector<std::pair<const DiracSpinor *, const DiracSpinor *>>
orbital_pairs(const std::vector<DiracSpinor> &a_orbs,
              const std::vector<DiracSpinor> &b_orbs,
              const DiracOperator::TensorOperator *h,
              const MEoptions &options) {

  std::vector<std::pair<const DiracSpinor *, const DiracSpinor *>> pairs;

  // The same list on both sides: each pair is calculated once. Two distinct
  // lists have no such pairing, so every pair is calculated
  const auto same_list = &a_orbs == &b_orbs;

  // Diagonal (only for even-parity operators)
  if (options.diagonal && (h->parity() == 1 || !same_list)) {
    for (const auto &a : a_orbs) {
      // Distinct lists: only the states that appear in both are 'diagonal'
      if (!same_list &&
          std::find(b_orbs.cbegin(), b_orbs.cend(), a) == b_orbs.cend())
        continue;
      if (h->isZero(a.kappa(), a.kappa()))
        continue;
      pairs.emplace_back(&a, &a);
    }
  }

  // Off-diagonal
  if (options.off_diagonal) {
    for (std::size_t ib = 0; ib < b_orbs.size(); ib++) {
      const auto &b = b_orbs.at(ib);
      for (std::size_t ia = 0; ia < a_orbs.size(); ia++) {
        const auto &a = a_orbs.at(ia);

        if (a == b)
          continue;
        if (h->isZero(a.kappa(), b.kappa()))
          continue;

        // Even-parity state on the right for odd-parity operators; upper
        // triangle for even-parity ones
        if (same_list && !options.calculate_both) {
          if (h->parity() == -1 && b.parity() == -1)
            continue;
          if (h->parity() == 1 && ib > ia)
            continue;
        }

        pairs.emplace_back(&a, &b);
      }
    }
  }

  return pairs;
}

//==============================================================================
std::vector<MEdata> matrix_elements(const std::vector<DiracSpinor> &a_orbs,
                                    const std::vector<DiracSpinor> &b_orbs,
                                    DiracOperator::TensorOperator *h,
                                    ExternalField::CorePolarisation *dV,
                                    MBPT::StructureRad *sr,
                                    const MEoptions &options,
                                    std::ostream &outstream) {

  const bool op_each = options.operator_omega == Frequency::transition;
  const bool rpa_each = options.rpa_omega == Frequency::transition;

  // The caller owns a fixed-frequency RPA. Solving it may not be wanted, so
  // this is a note, not an error
  if (!rpa_each && dV && dV->last_its() == 0) {
    fmt::print(outstream,
               "\nNote: RPA at a fixed frequency, but solve_core() "
               "has not been called (or broke):\n ** dV will be zero\n");
  }

  // With everything at a fixed frequency, the SR tables are built once here
  const bool freq_dep = h->freqDependantQ();
  if (sr && !((op_each && freq_dep) || rpa_each)) {
    sr->solve_core(h, dV);
  }

  std::vector<MEdata> mes;
  bool diagonal_set = false;
  for (const auto &[a, b] : orbital_pairs(a_orbs, b_orbs, h, options)) {

    const auto w_ab = a->en() - b->en();

    // The diagonal elements all share w = 0: set once
    if (a != b) {
      if (rpa_each && dV && options.print) {
        fmt::print(outstream,
                   "{} -> {}: w = {:.8f}\n RPA(w) : ", b->shortSymbol(),
                   a->shortSymbol(), w_ab);
      }
      solve_at_frequency(h, dV, sr, w_ab, options);
    } else if (!diagonal_set) {
      solve_at_frequency(h, dV, sr, 0.0, options);
      diagonal_set = true;
    }

    mes.push_back(matrix_element(*a, *b, h, dV, sr, options.type, w_ab,
                                 sr_frequency(dV, w_ab, options)));
  }

  return mes;
}

//==============================================================================
Coulomb::meTable<double> me_table(const std::vector<DiracSpinor> &a_orbs,
                                  const std::vector<DiracSpinor> &b_orbs,
                                  const DiracOperator::TensorOperator *h) {
  return me_table(a_orbs, b_orbs, h, nullptr, nullptr, 0.0);
}

//==============================================================================
Coulomb::meTable<double> me_table(const std::vector<DiracSpinor> &a_orbs,
                                  const std::vector<DiracSpinor> &b_orbs,
                                  const DiracOperator::TensorOperator *h,
                                  const ExternalField::CorePolarisation *dV,
                                  const MBPT::StructureRad *srn, double omega,
                                  int sr_n_max, bool sr_norm) {

  Coulomb::meTable<double> h_ab;

  const auto a_is_b = &a_orbs == &b_orbs;

  // Every element carries the RPA component at omega: the emission one
  // (conjugate) for a negative frequency
  const auto conj = omega < 0.0;

  // At omega = 0 (and the RPA solved there, or absent) the table is
  // symmetric: dV and the structure radiation are both Hermitian, so the
  // reversed pair is the symmetry conjugate. Otherwise it must be calculated:
  // dV(omega) of (b, a) is not the conjugate of that of (a, b), and the
  // structure radiation of (b, a) at omega is that of (a, b) at -omega
  const auto symmetric =
    omega == 0.0 && (dV == nullptr || dV->last_omega() == 0.0);

  // SR+N only between physical (low-n) states
  const auto sr_included = [&](const DiracSpinor &x, const DiracSpinor &y) {
    return srn && x.n() <= sr_n_max && y.n() <= sr_n_max;
  };

  const auto matrix_element = [&](const DiracSpinor &x, const DiracSpinor &y) {
    const auto tab = h->reducedME(x, y);
    const auto dv = dV ? dV->dV(x, y, conj) : 0.0;
    // Optionally without the norm
    const auto sr = !sr_included(x, y) ? 0.0 :
                    sr_norm            ? srn->srn(x, y, h, dV, omega) :
                                         srn->SR(x, y, omega);
    return tab + dv + sr;
  };

  // Each element to calculate is its own task, for an even load balance (the
  // cost is concentrated in the few SR+N elements). The reversed element is a
  // separate task, unless it follows by symmetry
  std::vector<std::pair<const DiracSpinor *, const DiracSpinor *>> elements;
  for (const auto &a : a_orbs) {
    for (const auto &b : b_orbs) {
      if (b < a && a_is_b)
        continue;
      if (h->isZero(a, b))
        continue;
      h_ab.add(a, b, 0.0);
      elements.emplace_back(&a, &b);
      if (a != b) {
        h_ab.add(b, a, 0.0);
        if (!symmetric) {
          elements.emplace_back(&b, &a);
        }
      }
    }
  }

  // SR+N elements cost far more than the rest: schedule them first
  std::stable_partition(elements.begin(), elements.end(),
                        [&](const auto &element) {
                          return sr_included(*element.first, *element.second);
                        });

#pragma omp parallel for schedule(dynamic)
  for (std::size_t i = 0; i < elements.size(); ++i) {
    const auto &x = *elements[i].first;
    const auto &y = *elements[i].second;

    const auto me = matrix_element(x, y);

    h_ab.update(x, y, me);
    if (symmetric && x != y) {
      h_ab.update(y, x, h->symm_sign(x, y) * me);
    }
  }
  return h_ab;
}

} // namespace Amplitudes
