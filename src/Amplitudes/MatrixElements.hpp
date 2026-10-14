#pragma once
#include "Coulomb/meTable.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include <iostream>
#include <string>
#include <utility>
#include <vector>
class DiracSpinor;
namespace ExternalField {
class CorePolarisation;
}
namespace MBPT {
class StructureRad;
}

//! Physical amplitudes and observables (matrix elements, second-order
//! amplitudes); testable functions, callable from any module.
namespace Amplitudes {

/*!
  @brief Result of a single matrix element calculation: the lowest-order
  value and each correction, separately.
  @details
  Holds \f$ \redmatel{a}{h}{b} \f$ and its corrections,
  \f[
    t^{\rm tot}_{ab} = t^{(0)}_{ab} + \delta V_{ab} + \delta t^{\rm SR}_{ab}
      + \delta t^{\rm Norm}_{ab}.
  \f]
  Everything is stored unscaled; the
  accessors apply the factor that converts the reduced matrix element to the
  requested form (reduced, stretched, or hyperfine constant); see
  DiracOperator::MatrixElementType.
*/
struct MEdata {
  //! State labels (shortSymbol); a is the bra: <a||h||b>
  std::string a{}, b{};
  //! Transition frequency, e_a - e_b (0 for diagonal)
  double omega{0.0};
  //! Factor converting reduced ME to requested MatrixElementType (1 for Reduced)
  double factor{1.0};
  //! Lowest-order reduced matrix element <a||h||b>
  double t0{0.0};
  //! RPA correction <a||dV||b> (0 if no RPA)
  double dv{0.0};
  //! True if an RPA correction was calculated (distinguishes from dv = 0)
  bool has_rpa{false};
  //! Structure radiation: top + bottom + centre diagrams (0 if no SR)
  double sr{0.0};
  //! Normalisation of states factor, StructureRad::norm_factor (0 if no SR)
  double norm_factor{0.0};
  //! Normalisation of states: norm_factor * (t0 + dv) (0 if no SR)
  double norm{0.0};

  //! Lowest-order value: factor * t0
  double value0() const { return factor * t0; }
  //! Lowest-order + RPA: factor * (t0 + dv)
  double value() const { return factor * (t0 + dv); }
  //! Full corrected value: factor * (t0 + dv + sr + norm)
  double total() const { return factor * (t0 + dv + sr + norm); }
};

/*!
  @brief Which frequency the operator, or the RPA, is evaluated at.
  @details
  There is no obvious default: the right choice depends on the calculation,
  so it must be given explicitly.

  - `transition`: the driver sets the frequency itself, to each pair's own
    transition frequency \f$ \omega_{ab} = \en_a - \en_b \f$; see
    @ref solve_at_frequency. This is the physically correct frequency
    for a transition.

  - `fixed`: the driver does not touch the frequency. The operator is assumed
    to already be at the intended frequency, and the RPA to already have been
    solved there, by the caller (the structure radiation tables are then
    built once, by the driver). Use
    for a fixed external field, or to reproduce a fixed-frequency
    calculation.
*/
enum class Frequency { transition, fixed };

/*!
  @brief Options for the matrix element drivers.
  @details
  The two frequency choices are independent and have no default: see
  @ref Frequency. The operator frequency usually matters most; the RPA (and
  structure radiation) frequency dependence is often very small, so `fixed`
  is usually adequate there.
*/
struct MEoptions {
  //! Frequency of the operator itself; see Frequency
  Frequency operator_omega;
  //! Frequency the RPA (and the structure radiation tables) are solved at;
  //! see Frequency. With transition and structure radiation, the SR tables
  //! are re-built at each transition frequency, which is slow
  Frequency rpa_omega;
  //! Calculate diagonal matrix elements (only for even-parity operators)
  bool diagonal{true};
  //! Calculate off-diagonal matrix elements
  bool off_diagonal{true};
  //! Calculate both <a||h||b> and <b||h||a>
  bool calculate_both{false};
  //! Form of matrix element: Reduced, Stretched, or HFConstant
  DiracOperator::MatrixElementType type{
    DiracOperator::MatrixElementType::Reduced};
  //! Maximum RPA iterations (1 corresponds to first-order RPA)
  int rpa_iterations{128};
  //! Print RPA solve progress
  bool print{true};

  //! The frequency choices must be made explicitly; there is no default
  MEoptions(Frequency t_operator_omega, Frequency t_rpa_omega)
    : operator_omega(t_operator_omega), rpa_omega(t_rpa_omega) {}
};

/*!
  @brief Sets the operator frequency, solves the RPA, and re-builds the
  structure radiation tables, for a transition at frequency omega, as the
  options ask. The one place the frequencies are handled.
  @details
  In this order:
  - with operator_omega = transition: the operator is set to
    \f$ +|\omega| \f$, since the RPA must be solved with the operator at
    positive frequency;
  - with rpa_omega = transition: the RPA is re-solved at @p omega (cleared
    first if the previous solution is poorly converged, or if first-order
    RPA was asked for, which must not iterate from a previous solution);
  - with operator_omega = transition: the operator is then set to the
    signed @p omega, its physical frequency for this transition;
  - with structure radiation, and either of the above: the SR tables are
    re-built, so that their internal lines carry the operator and the RPA
    as just set.

  With both frequencies fixed, nothing is touched: the caller has set the
  operator, solved the RPA, and built the SR tables.

  @param h       The operator (frequency updated if frequency-dependent).
  @param dV      RPA; nullptr for none.
  @param sr      Structure radiation; nullptr for none.
  @param omega   The transition frequency, signed: e_a - e_b for <a||h||b>.
  @param options See MEoptions (the two Frequency choices, rpa_iterations,
                 print).
*/
void solve_at_frequency(DiracOperator::TensorOperator *h,
                        ExternalField::CorePolarisation *dV,
                        MBPT::StructureRad *sr, double omega,
                        const MEoptions &options);

/*!
  @brief Single matrix element of h between states a and b, with its RPA and
  structure radiation corrections.
  @details
  Pure evaluation: assumes @p h is at the signed transition frequency, and
  @p dV and @p sr solved there (see @ref solve_at_frequency). The RPA
  correction is the component at the transition frequency: the emission one
  (\f$ \delta V^\dagger \f$) if @p omega is negative; see
  CorePolarisation::dV. The structure radiation denominators are taken at
  @p omega_sr: its frequency dependence is very weak, so when the RPA is at
  a fixed frequency this is usually the frequency the RPA was solved at
  (see @ref sr_frequency). The Brueckner orbital correction
  (MBPT::StructureRad::BO) is not included here.

  The MatrixElementType factor is calculated with @p h (it is purely
  angular), and stored in the returned MEdata rather than applied.

  @param a,b     States: <a||h||b>.
  @param h       The tensor operator.
  @param dV      RPA correction, already solved; nullptr for none.
  @param sr      Structure radiation, tables already built; nullptr for none.
  @param type    Form of matrix element (Reduced, Stretched, HFConstant).
  @param omega   Transition frequency, e_a - e_b: selects the RPA component,
                 and is recorded in the output.
  @param omega_sr  Frequency for the structure radiation denominators.
  @return MEdata holding t0, dv, sr, norm, factor, omega, labels.
*/
[[nodiscard]] MEdata
matrix_element(const DiracSpinor &a, const DiracSpinor &b,
               const DiracOperator::TensorOperator *h,
               const ExternalField::CorePolarisation *dV = nullptr,
               const MBPT::StructureRad *sr = nullptr,
               DiracOperator::MatrixElementType type =
                 DiracOperator::MatrixElementType::Reduced,
               double omega = 0.0, double omega_sr = 0.0);

/*!
  @brief The frequency for the structure radiation denominators of a
  transition at omega: the transition frequency if the RPA follows each
  transition, otherwise the frequency the RPA was solved at (zero if there
  is no RPA).
  @details
  The structure radiation follows the RPA: its frequency dependence is very
  weak, so it is taken at the same frequency.
*/
double sr_frequency(const ExternalField::CorePolarisation *dV, double omega,
                    const MEoptions &options);

/*!
  @brief Frequency-derivative term of <a||h||b> for a frequency-dependent
  operator: (dt/dw) dw2, the shift from the second-order energy shifts
  dw2 moving the transition frequency.
  @details
  \f$ \delta t = ({\rm d}t/{\rm d}\omega)\,\delta\omega^{(2)} \f$.
  Off-diagonal elements only; zero for a frequency-independent operator.
  The derivative is taken numerically, which is why the operator is not
  const. It is left as it was, at @p omega.

  @param a,b     States: <a||h||b>.
  @param h       The operator, at the signed transition frequency @p omega.
  @param sr      Structure radiation (holds Sigma for the energy shifts).
  @param omega   Transition frequency, e_a - e_b.
*/
[[nodiscard]] double frequency_derivative_term(const DiracSpinor &a,
                                               const DiracSpinor &b,
                                               DiracOperator::TensorOperator *h,
                                               const MBPT::StructureRad *sr,
                                               double omega);

/*!
  @brief The pairs of orbitals (a, b) of <a||h||b> that the drivers
  calculate, in order.
  @details
  Diagonal elements first (only for even-parity operators, and only for
  states in both lists), then off-diagonal; pairs forbidden by the selection
  rules are skipped. When @p a_orbs and @p b_orbs are the same list, each
  pair is calculated once: for odd-parity operators only elements with the
  even-parity state on the right, for even-parity operators only the upper
  triangle (unless calculate_both). Two distinct lists have no such pairing,
  so every pair is included.

  @param a_orbs  Bra states (a).
  @param b_orbs  Ket states (b).
  @param h       The tensor operator (selection rules).
  @param options diagonal, off_diagonal, calculate_both.
  @return Pointers into @p a_orbs and @p b_orbs: {a, b} for each element.
*/
[[nodiscard]] std::vector<std::pair<const DiracSpinor *, const DiracSpinor *>>
orbital_pairs(const std::vector<DiracSpinor> &a_orbs,
              const std::vector<DiracSpinor> &b_orbs,
              const DiracOperator::TensorOperator *h, const MEoptions &options);

/*!
  @brief Matrix elements of h for all allowed pairs from two lists of
  orbitals, with optional RPA and structure radiation; owns all frequency
  updates and RPA solves.
  @details
  For each pair of @ref orbital_pairs: @ref solve_at_frequency, then
  @ref matrix_element. Nothing is printed beyond the RPA progress (if
  options.print, which also prints each transition as the RPA is
  re-solved).

  @param a_orbs    Bra states (index a).
  @param b_orbs    Ket states (index b).
  @param h         The tensor operator.
  @param dV        RPA. nullptr for no RPA.
  @param sr        Structure radiation (mutated: its tables are built).
                   nullptr for none.
  @param options   See MEoptions.
  @param outstream Stream for progress output.
  @return Vector of MEdata, one per calculated matrix element.
*/
[[nodiscard]] std::vector<MEdata>
matrix_elements(const std::vector<DiracSpinor> &a_orbs,
                const std::vector<DiracSpinor> &b_orbs,
                DiracOperator::TensorOperator *h,
                ExternalField::CorePolarisation *dV, MBPT::StructureRad *sr,
                const MEoptions &options, std::ostream &outstream = std::cout);

/*!
  @brief Matrix elements of h for all pairs from a single list of orbitals.
  @details
  Convenience overload; calls matrix_elements(orbs, orbs, ...) with both
  bra and ket taken from @p orbs, so each pair is calculated once.
*/
[[nodiscard]] inline std::vector<MEdata>
matrix_elements(const std::vector<DiracSpinor> &orbs,
                DiracOperator::TensorOperator *h,
                ExternalField::CorePolarisation *dV, MBPT::StructureRad *sr,
                const MEoptions &options, std::ostream &outstream = std::cout) {
  return matrix_elements(orbs, orbs, h, dV, sr, options, outstream);
}

//==============================================================================

/*!
  @brief Builds a lookup table of reduced matrix elements <a||h||b>.
  @details
  Fills and returns a `Coulomb::meTable<double>` with the bare reduced
  matrix elements \f$ \redmatel{a}{h}{b} \f$ for all non-zero pairs from
  @p a_orbs and @p b_orbs. The symmetry-conjugate \f$ \redmatel{b}{h}{a} \f$
  is also stored, via `symm_sign()`. Filled with OpenMP parallelisation.
  @p h must already be at the intended frequency. For a table with RPA (and
  structure radiation), use the overload taking @p omega.

  @param a_orbs  Bra states.
  @param b_orbs  Ket states.
  @param h       Pointer to the (const) tensor operator.

  @return meTable containing t_ab for all non-zero pairs (and conjugates).
*/
[[nodiscard]] Coulomb::meTable<double>
me_table(const std::vector<DiracSpinor> &a_orbs,
         const std::vector<DiracSpinor> &b_orbs,
         const DiracOperator::TensorOperator *h);

/*!
  @brief Builds a matrix element table for a single set of orbitals.
  @details
  Convenience overload; calls me_table(a_orbs, a_orbs, ...) with both
  bra and ket taken from @p a_orbs.
*/
[[nodiscard]] inline Coulomb::meTable<double>
me_table(const std::vector<DiracSpinor> &a_orbs,
         const DiracOperator::TensorOperator *h) {
  return me_table(a_orbs, a_orbs, h);
}

/*!
  @brief Builds a table of reduced matrix elements, including structure
  radiation and (optionally) the normalisation of states.
  @details
  As above, but each element also carries the second-order corrections,
  \f[ t_{ab} = \redmatel{a}{h}{b} + \delta V_{ab} + \delta_{\rm SR}^{ab}, \f]
  and every element is at the one frequency @p omega: the RPA correction is
  \f$ \delta V(\omega) \f$ throughout (the emission component,
  \f$ \delta V^\dagger \f$, for negative @p omega; see
  CorePolarisation::dV), and the structure radiation is at @p omega,
  absorbed by the ket. The reversed pair \f$ (b, a) \f$ is then
  not the symmetry conjugate of \f$ (a, b) \f$, and is calculated as well
  (except at \f$ \omega = 0 \f$, where it is).

  @param a_orbs  Bra states.
  @param b_orbs  Ket states.
  @param h       Pointer to the (const) tensor operator.
  @param dV      Optional RPA correction. If nullptr, not applied.
  @param srn     Structure radiation/normalisation. If nullptr, not applied
                 (and the table is as the plain overload above).
  @param omega   Frequency of the table: that of the operator component the
                 table stands for. Sets the sign of the RPA correction, and
                 the structure radiation denominators.
  @param sr_n_max  SR+N is applied only to pairs with both n <= @p sr_n_max.
                 SR+N is meaningful only between physical states, so this
                 limits it to the low-n part of a large basis, where the
                 states are not cavity states. Does not affect the internal
                 lines of the diagrams (see MBPT::StructureRad) [999].
  @param sr_norm If false, only the structure radiation is added, not the
                 normalisation of states [true].

  @return meTable containing t_ab for all non-zero pairs, both orderings.
*/
[[nodiscard]] Coulomb::meTable<double>
me_table(const std::vector<DiracSpinor> &a_orbs,
         const std::vector<DiracSpinor> &b_orbs,
         const DiracOperator::TensorOperator *h,
         const ExternalField::CorePolarisation *dV,
         const MBPT::StructureRad *srn, double omega, int sr_n_max = 999,
         bool sr_norm = true);

/*!
  @brief Builds a SR+N matrix element table for a single set of orbitals.
  @details
  Convenience overload; calls me_table(a_orbs, a_orbs, ...) with both
  bra and ket taken from @p a_orbs.
*/
[[nodiscard]] inline Coulomb::meTable<double>
me_table(const std::vector<DiracSpinor> &a_orbs,
         const DiracOperator::TensorOperator *h,
         const ExternalField::CorePolarisation *dV,
         const MBPT::StructureRad *srn, double omega, int sr_n_max = 999,
         bool sr_norm = true) {
  return me_table(a_orbs, a_orbs, h, dV, srn, omega, sr_n_max, sr_norm);
}

//==============================================================================

} // namespace Amplitudes
