#include "Amplitudes/MatrixElements.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "IO/ChronoTimer.hpp"
#include "MBPT/StructureRad.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/ostream.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace Amplitudes {

//==============================================================================
MEdata matrix_element(const DiracSpinor &a, const DiracSpinor &b,
                      const DiracOperator::TensorOperator *h,
                      const DiracOperator::TensorOperator *h_minus,
                      const ExternalField::CorePolarisation *dV,
                      DiracOperator::MatrixElementType type, double omega) {

  // t_- (at negative frequency) for sign-sensitive freq-dependent operators
  const auto *h_pm = (omega < 0.0 && h_minus) ? h_minus : h;

  MEdata me;
  me.a = a.shortSymbol();
  me.b = b.shortSymbol();
  me.omega = omega;
  me.factor = h->matel_factor(type, a, b);
  me.t0 = h_pm->reducedME(a, b);
  me.dv = dV ? dV->dV(a, b) : 0.0;
  me.has_rpa = (dV != nullptr);
  return me;
}

//==============================================================================
void set_operator_frequency(DiracOperator::TensorOperator *h,
                            DiracOperator::TensorOperator *h_minus,
                            double omega) {
  if (h == nullptr || !h->freqDependantQ())
    return;
  h->updateFrequency(std::abs(omega));
  if (h_minus) {
    h_minus->updateFrequency(-std::abs(omega));
  }
}

//==============================================================================
std::vector<MEdata> matrix_elements(const std::vector<DiracSpinor> &a_orbs,
                                    const std::vector<DiracSpinor> &b_orbs,
                                    DiracOperator::TensorOperator *h,
                                    DiracOperator::TensorOperator *h_minus,
                                    ExternalField::CorePolarisation *dV,
                                    const MEoptions &options,
                                    std::ostream &outstream) {

  std::vector<MEdata> mes;

  // The same list on both sides: each pair is calculated once. Two distinct
  // lists have no such pairing, so every pair is calculated
  const auto same_list = &a_orbs == &b_orbs;

  // The operator, and the RPA, are each either updated here at every
  // transition frequency, or left exactly as the caller set them
  const bool op_each = options.operator_omega == Frequency::transition;
  const bool rpa_each = options.rpa_omega == Frequency::transition;

  // The caller owns a fixed-frequency RPA. Solving it may not be wanted, so
  // this is a note, not an error
  if (!rpa_each && dV && dV->last_its() == 0) {
    fmt::print(outstream, "\nNote: RPA at Frequency::fixed, and solve_core() "
                          "has not been called: dV will be zero\n");
  }

  // Re-solve RPA at given frequency; clear first if previous solution poorly
  // converged, or if first-order RPA (must not iterate from prior solution)
  const auto resolve_rpa = [&](double w) {
    if (!dV)
      return;
    if (dV->last_eps() > 1.0e-5 || options.rpa_iterations == 1 ||
        std::isnan(dV->last_eps())) {
      dV->clear();
    }
    dV->solve_core(w, options.rpa_iterations, options.print);
  };

  //----------------------------------------------------
  // Diagonal (only for even-parity operators): w = 0
  if (options.diagonal && (h->parity() == 1 || !same_list)) {

    if (op_each) {
      set_operator_frequency(h, h_minus, 0.0);
    }
    if (rpa_each) {
      resolve_rpa(0.0);
    }

    for (const auto &a : a_orbs) {
      // Distinct lists: only the states that appear in both are 'diagonal'
      if (!same_list &&
          std::find(b_orbs.cbegin(), b_orbs.cend(), a) == b_orbs.cend())
        continue;
      if (h->isZero(a.kappa(), a.kappa()))
        continue;
      mes.push_back(matrix_element(a, a, h, h_minus, dV, options.type, 0.0));
    }
  }

  //----------------------------------------------------
  // Off-diagonal:
  if (options.off_diagonal) {
    for (std::size_t ib = 0; ib < b_orbs.size(); ib++) {
      const auto &b = b_orbs.at(ib);
      for (std::size_t ia = 0; ia < a_orbs.size(); ia++) {
        const auto &a = a_orbs.at(ia);

        if (a == b)
          continue;
        if (h->isZero(a.kappa(), b.kappa()))
          continue;

        // Ensure even-parity state on right for odd-parity operators.
        // Only for a single list: with two lists, every pair is calculated
        if (same_list) {
          if (h->parity() == -1) {
            if (!options.calculate_both && b.parity() == -1)
              continue;
          } else {
            if (!options.calculate_both && ib > ia)
              continue;
          }
        }

        const auto w_ab = a.en() - b.en();

        if (op_each) {
          set_operator_frequency(h, h_minus, w_ab);
        }
        if (rpa_each && dV) {
          if (options.print) {
            fmt::print(outstream,
                       "{} -> {}: w = {:.8f}\n RPA(w) : ", b.shortSymbol(),
                       a.shortSymbol(), w_ab);
          }
          resolve_rpa(w_ab);
        }

        mes.push_back(matrix_element(a, b, h, h_minus, dV, options.type, w_ab));
      }
    }
  }

  return mes;
}

//==============================================================================
Coulomb::meTable<double> me_table(const std::vector<DiracSpinor> &a_orbs,
                                  const std::vector<DiracSpinor> &b_orbs,
                                  const DiracOperator::TensorOperator *h,
                                  const ExternalField::CorePolarisation *dV) {
  return me_table(a_orbs, b_orbs, h, dV, nullptr, 0.0);
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

  for (const auto &a : a_orbs) {
    for (const auto &b : b_orbs) {
      if (b < a && a_is_b)
        continue;
      if (h->isZero(a, b))
        continue;
      h_ab.add(a, b, 0.0);
      if (a != b) {
        h_ab.add(b, a, 0.0);
      }
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (std::size_t i = 0; i < a_orbs.size(); ++i) {
    const auto &a = a_orbs[i];
    for (const auto &b : b_orbs) {

      if (b < a && a_is_b)
        continue;

      if (h->isZero(a, b))
        continue;

      const auto tab = h->reducedME(a, b);
      const auto dv = dV ? dV->dV(a, b) : 0.0;
      // SR+N only between physical (low-n) states; optionally without the norm
      const auto do_sr = srn && a.n() <= sr_n_max && b.n() <= sr_n_max;
      const auto sr = !do_sr  ? 0.0 :
                      sr_norm ? srn->srn(a, b, h, dV, omega) :
                                srn->SR(a, b, omega);

      const auto me = tab + dv + sr;

      h_ab.update(a, b, me);
      if (a != b) {
        h_ab.update(b, a, h->symm_sign(a, b) * me);
      }
    }
  }
  return h_ab;
}


//==============================================================================
//==============================================================================
Coulomb::meTable<double> me_table_Jpnc(const std::vector<DiracSpinor> &a_orbs,
  const std::vector<DiracSpinor> &b_orbs,
  const DiracOperator::TensorOperator *h,
  const ExternalField::CorePolarisation *dV) {

Coulomb::meTable<double> h_ab;

const auto a_is_b = &a_orbs == &b_orbs;

for (const auto &a : a_orbs) {
for (const auto &b : b_orbs) {
if (b < a && a_is_b)
continue;
if (h->isZero(a, b))
continue;
h_ab.add(a, b, 0.0);
if (a != b) {
h_ab.add(b, a, 0.0);
}
}
}

#pragma omp parallel for schedule(dynamic)
for (std::size_t i = 0; i < a_orbs.size(); ++i) {
const auto &a = a_orbs[i];
for (const auto &b : b_orbs) {

if (b < a && a_is_b)
continue;

if (h->isZero(a, b))
continue;


//Johnson factor 
double tja = a.twoj();
const auto factor = 1.0/sqrt(tja+1);
const auto tab = h->reducedME(a, b);
const auto tab_J = tab*factor;

const auto dv = dV ? dV->dV(a, b) : 0.0;


const auto me = tab_J + dv;

h_ab.update(a, b, me);
if (a != b) {
h_ab.update(b, a, h->symm_sign(a, b) * me);
}
}
}
return h_ab;
}



//==============================================================================
std::vector<SRNdata> sr_matrix_elements(const std::vector<DiracSpinor> &orbs,
                                        DiracOperator::TensorOperator *h,
                                        MBPT::StructureRad *sr,
                                        ExternalField::CorePolarisation *dV,
                                        const SRNoptions &options,
                                        std::ostream &outstream) {
  assert(sr != nullptr);

  std::vector<SRNdata> out;

  const bool freq_dep = h->freqDependantQ();
  const bool op_each = options.operator_omega == Frequency::transition;
  const bool rpa_each = options.rpa_omega == Frequency::transition;

  // The caller owns a fixed-frequency RPA. Solving it may not be wanted, so
  // this is a note, not an error
  if (!rpa_each && dV && dV->last_its() == 0) {
    fmt::print(outstream, "\nNote: RPA at Frequency::fixed, and solve_core() "
                          "has not been called: dV will be zero\n");
  }

  // The structure radiation follows the RPA: its frequency dependence is
  // very weak, so when the RPA is fixed, the SR is taken at the frequency
  // the RPA was solved at (zero if there is no RPA)
  const auto fixed_omega = dV ? dV->last_omega() : 0.0;

  // The SR matrix element tables hold <a||h+dV||b> for the internal lines,
  // so they must be built after the operator and the RPA are at the right
  // frequency; this is done inside the loop below. The frequency they were
  // last built at, so they are re-built only when it changes:
  std::optional<double> tables_at{};

  // Diagonal first, then off-diagonal
  for (const auto diag : {true, false}) {

    if (diag && !options.diagonal)
      continue;
    if (!diag && !options.off_diagonal)
      continue;

    for (std::size_t ib = 0; ib < orbs.size(); ib++) {
      const auto &v = orbs.at(ib);
      for (std::size_t ia = 0; ia < orbs.size(); ia++) {
        const auto &w = orbs.at(ia);

        if (h->isZero(w.kappa(), v.kappa()))
          continue;
        // Ensure even-parity state on right for odd-parity operators
        if (h->parity() == -1) {
          if (!options.calculate_both && v.parity() == -1)
            continue;
        } else {
          if (!options.calculate_both && ib > ia)
            continue;
        }
        if (diag != (v == w))
          continue;

        std::optional<IO::ChronoTimer> timer;
        if (options.print) {
          timer.emplace("time");
        }

        // Transition frequency (operator), and the RPA/SR frequency
        const auto w_ab = w.en() - v.en();
        const auto ww_sr = rpa_each ? w_ab : fixed_omega;
        if (options.print) {
          fmt::print(outstream, "\n{} -> {}: w = {:.8f}\n", v.shortSymbol(),
                     w.shortSymbol(), w_ab);
        }

        const auto factor = h->matel_factor(options.type, w, v);

        // The operator itself is at the (physical) transition frequency
        if (op_each && freq_dep) {
          h->updateFrequency(w_ab);
        }
        // Re-solve the RPA at this transition frequency
        if (rpa_each && dV) {
          if (dV->last_eps() > 1.0e-3 || std::isnan(dV->last_eps())) {
            dV->clear();
          }
          dV->solve_core(w_ab, 100, options.print);
        }
        // Build the SR tables, now that the operator and the RPA are set for
        // this pair. Only re-built when their frequency has moved
        const auto table_w =
          (op_each && freq_dep) || rpa_each ? w_ab : fixed_omega;
        if (!tables_at || *tables_at != table_w) {
          sr->solve_core(h, dV);
          tables_at = table_w;
        }

        SRNdata me;
        me.a = w.shortSymbol();
        me.b = v.shortSymbol();
        me.omega = w_ab;
        me.factor = factor;
        me.t0 = h->reducedME(w, v);
        me.dv = dV ? dV->dV(w, v) : 0.0;
        me.has_rpa = (dV != nullptr);

        if (options.print) {
          fmt::print(outstream, "{:8s}  {:12.5e}\n", "t0", factor * me.t0);
          if (dV) {
            fmt::print(outstream, "{:8s}  {:12.5e}\n", "dV", factor * me.dv);
            fmt::print(outstream, "{:8s}  {:12.5e}\n", "RPA", me.value0());
          }
          outstream << std::flush;
        }

        // "Top" + "Bottom" + "Centre" SR diagrams:
        me.sr = sr->SR(w, v, ww_sr);
        if (options.print) {
          fmt::print(outstream, "{:8s}  {:12.5e}\n", "SR", factor * me.sr);
          outstream << std::flush;
        }

        // Normalisation of states:
        const auto f_norm =
          sr->f_norm(w) + (w == v ? sr->f_norm(w) : sr->f_norm(v));
        me.norm = f_norm * (me.t0 + me.dv);
        if (options.print) {
          fmt::print(outstream, "{:8s}  {:12.5e}\n", "Norm(0)",
                     factor * me.norm);
          fmt::print(outstream, "{:8s}  {:12.5e}\n", "Norm(SR)",
                     factor * f_norm * (me.t0 + me.dv + me.sr));
          outstream << std::flush;
        }

        // Brueckner orbital correction (only when the legs are not already
        // Brueckner orbitals):
        if (options.include_bo) {
          me.bo = sr->BO(w, v);
          if (options.print) {
            fmt::print(outstream, "{:8s}  {:12.5e}\n", "BO", factor * me.bo);
          }

          // Frequency-derivative term: the second-order energy shifts move
          // the transition frequency; for a freq-dependent operator this
          // shifts the matrix element: dT = (dt/dw) * dw^(2)
          if (freq_dep && op_each && w != v) {
            const auto dw_2 = sr->Sigma_vw(w, w) - sr->Sigma_vw(v, v);
            const auto del = 0.05 * std::abs(w_ab);
            h->updateFrequency(w_ab + del);
            const auto h_plus = h->reducedME(w, v);
            h->updateFrequency(w_ab - del);
            const auto h_minus = h->reducedME(w, v);
            h->updateFrequency(w_ab);
            const auto T_deriv = (h_plus - h_minus) / (2 * del) * dw_2;
            if (options.print) {
              fmt::print(outstream, "+ {:7s} {:12.5e}\n", "Tderiv",
                         factor * T_deriv);
            }
            me.bo += T_deriv;
          }
        }

        if (options.print) {
          fmt::print(outstream, "{:8s}  {:12.5e}\n", "Total", me.total());
          outstream << std::flush;
        }

        out.push_back(me);
      }
    }
  }

  return out;
}

} // namespace Amplitudes
