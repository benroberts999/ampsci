#pragma once
#include "Coulomb/meTable.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include <iostream>
#include <optional>
#include <string>
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
  @brief Result of a single matrix element calculation.
  @details
  Holds \f$ \redmatel{a}{h}{b} \f$: the lowest-order value, the RPA (core
  polarisation) correction, the transition frequency, and the factor that
  converts the reduced matrix element to the requested form (reduced,
  stretched, or hyperfine constant); see DiracOperator::MatrixElementType.

  Everything is stored unscaled; value() and value0() apply the factor.
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

  //! Full value: factor * (t0 + dv)
  double value() const { return factor * (t0 + dv); }
  //! Lowest-order value: factor * t0
  double value0() const { return factor * t0; }
};

/*!
  @brief Which frequency the operator, or the RPA, is evaluated at.
  @details
  There is no obvious default: the right choice depends on the calculation,
  so it must be given explicitly.

  - `transition`: the driver sets the frequency itself, to each pair's own
    transition frequency \f$ \omega_{ab} = \en_a - \en_b \f$. The operator is
    updated (h at \f$ +|\omega_{ab}| \f$, h_minus at \f$ -|\omega_{ab}| \f$),
    and the RPA is re-solved, at every element. This is the physically
    correct frequency for a transition.

  - `fixed`: the driver does not touch the frequency. The operator is assumed
    to already be at the intended frequency, and the RPA to already have been
    solved there, by the caller. Use for a fixed external field, or to
    reproduce a fixed-frequency calculation.
*/
enum class Frequency { transition, fixed };

/*!
  @brief Options for the matrix_elements() list driver.
  @details
  The two frequency choices are independent and have no default: see
  @ref Frequency. The operator frequency usually matters most; the RPA
  frequency dependence is often very small, so `fixed` is usually adequate
  there.
*/
struct MEoptions {
  //! Frequency of the operator itself; see Frequency
  Frequency operator_omega;
  //! Frequency the RPA is solved at; see Frequency
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
  @brief Sets the \f$ t_\pm \f$ operator pair to the frequency w.
  @details
  Sets @p h to \f$ +|\omega| \f$ and, if given, @p h_minus to
  \f$ -|\omega| \f$; does nothing for frequency-independent operators, for
  which updateFrequency() must not be called.
  Call this before @ref matrix_elements when using Frequency::fixed.

  @note Only for the \f$ t_\pm \f$ pair. Where there is no \f$ t_- \f$
        (e.g. @ref sr_matrix_elements), the operator takes the signed
        frequency directly, not \f$ |\omega| \f$.
*/
void set_operator_frequency(DiracOperator::TensorOperator *h,
                            DiracOperator::TensorOperator *h_minus,
                            double omega);

/*!
  @brief Single matrix element of h between states a and b, with optional RPA.
  @details
  Pure evaluation: assumes @p h, @p h_minus, and @p dV are already at the
  correct frequency. If @p omega is negative and @p h_minus is given,
  @p h_minus is used for the matrix element (sign-sensitive
  frequency-dependent operators, e.g. E1v: h holds \f$ t_+ \f$ at
  \f$ +|\omega| \f$, h_minus holds \f$ t_- \f$ at \f$ -|\omega| \f$).

  The MatrixElementType factor is calculated with @p h (it is purely
  angular), and stored in the returned MEdata rather than applied.

  @param a,b     States: <a||h||b>.
  @param h       The tensor operator.
  @param h_minus Operator at negative frequency; nullptr if not required.
  @param dV      RPA correction, already solved; nullptr for none.
  @param type    Form of matrix element (Reduced, Stretched, HFConstant).
  @param omega   Transition frequency (only selects h vs h_minus, and is
                 recorded in the output).
  @return MEdata holding t0, dv, factor, omega, and labels.
*/
[[nodiscard]] MEdata
matrix_element(const DiracSpinor &a, const DiracSpinor &b,
               const DiracOperator::TensorOperator *h,
               const DiracOperator::TensorOperator *h_minus = nullptr,
               const ExternalField::CorePolarisation *dV = nullptr,
               DiracOperator::MatrixElementType type =
                 DiracOperator::MatrixElementType::Reduced,
               double omega = 0.0);

/*!
  @brief Matrix elements of h for all allowed pairs from two lists of
  orbitals, with optional RPA; owns all frequency updates and RPA solves.
  @details
  Calculates \f$ \redmatel{a}{h}{b} \f$ for each pair allowed by the
  selection rules, with the bra states taken from @p a_orbs and the ket
  states from @p b_orbs, diagonal first, then off-diagonal.

  Selection rules: pairs with isZero() are skipped; diagonal elements only
  for even-parity operators.

  When @p a_orbs and @p b_orbs are the same list (as in the single-list
  overload below), each pair is calculated once: for odd-parity operators,
  only elements with the even-parity state on the right are included
  (unless calculate_both); for even-parity operators, only the upper
  triangle. Two distinct lists have no such pairing, so every pair is
  calculated and calculate_both has no effect.

  Frequency handling (see @ref Frequency): with `transition`, the operator
  is updated at each pair's transition frequency, @p h at
  \f$ +|\omega_{ab}| \f$ and @p h_minus at \f$ -|\omega_{ab}| \f$, and the
  RPA is re-solved there (cleared first when poorly converged, or when
  rpa_iterations is 1 so that first-order RPA is not iterated from a
  previous solution). With `fixed`, neither is touched: the caller must set
  the operator frequency (see @ref set_operator_frequency) and solve the
  RPA before calling.

  @param a_orbs    Bra states (index a).
  @param b_orbs    Ket states (index b).
  @param h         The tensor operator.
  @param h_minus   Operator at negative frequency (e.g. a clone of @p h for
                   E1v); nullptr if not required.
  @param dV        RPA. nullptr for no RPA.
  @param options   See MEoptions.
  @param outstream Stream for progress output.
  @return Vector of MEdata, one per calculated matrix element.
*/
[[nodiscard]] std::vector<MEdata> matrix_elements(
  const std::vector<DiracSpinor> &a_orbs,
  const std::vector<DiracSpinor> &b_orbs, DiracOperator::TensorOperator *h,
  DiracOperator::TensorOperator *h_minus, ExternalField::CorePolarisation *dV,
  const MEoptions &options, std::ostream &outstream = std::cout);

/*!
  @brief Matrix elements of h for all pairs from a single list of orbitals.
  @details
  Convenience overload; calls matrix_elements(orbs, orbs, ...) with both
  bra and ket taken from @p orbs, so each pair is calculated once.
*/
[[nodiscard]] inline std::vector<MEdata> matrix_elements(
  const std::vector<DiracSpinor> &orbs, DiracOperator::TensorOperator *h,
  DiracOperator::TensorOperator *h_minus, ExternalField::CorePolarisation *dV,
  const MEoptions &options, std::ostream &outstream = std::cout) {
  return matrix_elements(orbs, orbs, h, h_minus, dV, options, outstream);
}

//==============================================================================

/*!
  @brief Builds a lookup table of reduced matrix elements <a||h||b>.
  @details
  Fills and returns a `Coulomb::meTable<double>` with reduced matrix elements
  \f[ t_{ab} = \redmatel{a}{h}{b} + \delta V_{ab} \f]
  for all non-zero pairs from @p a_orbs and @p b_orbs.

  The symmetry-conjugate \f$ \redmatel{b}{h}{a} \f$ is also stored, via
  `symm_sign()`. Filled with OpenMP parallelisation.

  This is a pure table builder: @p h must already be at the intended
  frequency, and @p dV already solved there.

  @param a_orbs  Bra states.
  @param b_orbs  Ket states.
  @param h       Pointer to the (const) tensor operator.
  @param dV      Optional RPA correction. If nullptr, not applied.

  @return meTable containing t_ab for all non-zero pairs (and conjugates).
*/
[[nodiscard]] Coulomb::meTable<double>
me_table(const std::vector<DiracSpinor> &a_orbs,
         const std::vector<DiracSpinor> &b_orbs,
         const DiracOperator::TensorOperator *h,
         const ExternalField::CorePolarisation *dV = nullptr);

/*!
  @brief Builds a matrix element table for a single set of orbitals.
  @details
  Convenience overload; calls me_table(a_orbs, a_orbs, ...) with both
  bra and ket taken from @p a_orbs.
*/
[[nodiscard]] inline Coulomb::meTable<double>
me_table(const std::vector<DiracSpinor> &a_orbs,
         const DiracOperator::TensorOperator *h,
         const ExternalField::CorePolarisation *dV = nullptr) {
  return me_table(a_orbs, a_orbs, h, dV);
}

/*!
  @brief Builds a table of reduced matrix elements, including structure
  radiation and (optionally) the normalisation of states.
  @details
  As above, but each element also carries the second-order corrections,
  \f[ t_{ab} = \redmatel{a}{h}{b} + \delta V_{ab} + \delta_{\rm SR}^{ab}. \f]

  @param a_orbs  Bra states.
  @param b_orbs  Ket states.
  @param h       Pointer to the (const) tensor operator.
  @param dV      Optional RPA correction. If nullptr, not applied.
  @param srn     Structure radiation/normalisation. If nullptr, not applied
                 (and the table is as the plain overload above).
  @param omega   Frequency for the structure radiation denominators. Their
                 frequency dependence is very weak, so in practice this is
                 usually taken as the frequency the RPA was solved at.
  @param sr_n_max  SR+N is applied only to pairs with both n <= @p sr_n_max.
                 SR+N is meaningful only between physical states, so this
                 limits it to the low-n part of a large basis, where the
                 states are not cavity states. Does not affect the internal
                 lines of the diagrams (see MBPT::StructureRad) [999].
  @param sr_norm If false, only the structure radiation is added, not the
                 normalisation of states [true].

  @return meTable containing t_ab for all non-zero pairs (and conjugates).
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

Coulomb::meTable<double> me_table_Jpnc(const std::vector<DiracSpinor> &a_orbs,
  const std::vector<DiracSpinor> &b_orbs,
  const DiracOperator::TensorOperator *h,
  const ExternalField::CorePolarisation *dV=nullptr);


//==============================================================================

/*!
  @brief Result of a matrix element calculation with second-order MBPT
  corrections: structure radiation, normalisation, Brueckner orbital.
  @details
  Everything is stored unscaled; the accessors apply the MatrixElementType
  factor. The total corrected matrix element is
  \f[
    t^{\rm tot}_{ab} = t^{(0)}_{ab} + \delta V_{ab} + \delta t^{\rm SR}_{ab}
      + \delta t^{\rm Norm}_{ab} + \delta t^{\rm BO}_{ab}.
  \f]
*/
struct SRNdata {
  //! State labels (shortSymbol); a is the bra: <a||h||b>
  std::string a{}, b{};
  //! Frequency the corrections were evaluated at
  double omega{0.0};
  //! Factor converting reduced ME to requested MatrixElementType
  double factor{1.0};
  //! Lowest-order reduced matrix element <a||h||b>
  double t0{0.0};
  //! RPA correction (0 if no RPA)
  double dv{0.0};
  //! Structure radiation (top + bottom + centre diagrams)
  double sr{0.0};
  //! Normalisation of states: (f_norm_a + f_norm_b) * (t0 + dv)
  double norm{0.0};
  //! Brueckner orbital correction (0 if legs are already Brueckner);
  //! includes the frequency-derivative term for freq-dependent operators
  double bo{0.0};
  //! True if an RPA correction was calculated
  bool has_rpa{false};

  //! Lowest-order + RPA: factor * (t0 + dv)
  double value0() const { return factor * (t0 + dv); }
  //! Full corrected value: factor * (t0 + dv + sr + norm + bo)
  double total() const { return factor * (t0 + dv + sr + norm + bo); }
};

/*!
  @brief Options for sr_matrix_elements()
  @details
  As MEoptions; see @ref Frequency. The structure radiation follows the RPA
  choice: its frequency dependence is very weak, so with `fixed` the SR
  denominators are evaluated at the frequency the RPA was solved at
  (zero if there is no RPA).
*/
struct SRNoptions {
  //! Frequency of the operator itself; see Frequency
  Frequency operator_omega;
  //! Frequency the RPA, and hence the structure radiation, is evaluated at.
  //! With transition, the SR matrix element tables are re-built at each
  //! transition frequency, which is slow
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
  //! Include the Brueckner orbital correction (set false when the external
  //! legs are already Brueckner orbitals)
  bool include_bo{true};
  //! Print per-element progress (SR is slow)
  bool print{true};

  //! The frequency choices must be made explicitly; there is no default
  SRNoptions(Frequency t_operator_omega, Frequency t_rpa_omega)
    : operator_omega(t_operator_omega), rpa_omega(t_rpa_omega) {}
};

/*!
  @brief Matrix elements with structure radiation, normalisation, and
  Brueckner orbital corrections, for all pairs from a list of orbitals.
  @details
  For each pair allowed by the selection rules (as @ref matrix_elements),
  evaluates the lowest-order matrix element, RPA, and the second-order MBPT
  corrections via MBPT::StructureRad: SR (top+bottom+centre diagrams),
  normalisation of states, and (optionally) the Brueckner orbital
  correction.

  Frequency handling (see @ref Frequency): with operator_omega = transition,
  the operator is evaluated at each pair's transition frequency, where its
  frequency dependence is important, and its BO term then includes the
  frequency-derivative correction
  \f$ ({\rm d}t/{\rm d}\omega)\,\delta\omega^{(2)} \f$ for off-diagonal
  elements. With rpa_omega = transition, the RPA and the SR tables are
  re-solved at each transition frequency; with fixed, the caller must have
  solved the RPA already, and the SR denominators are taken at the frequency
  it was solved at (zero if there is no RPA).

  The caller constructs (and owns) the MBPT::StructureRad object, which
  holds the basis, Qk integrals, and screening options; solve_core is called
  here, and re-called whenever the frequency changes.

  @param orbs      Orbitals for the external legs; all pairs considered.
  @param h         The tensor operator.
  @param sr        StructureRad object (mutated: solve_core is called).
  @param dV        RPA. nullptr for no RPA.
  @param options   See SRNoptions.
  @param outstream Stream for per-element progress output.
  @return Vector of SRNdata, one per calculated matrix element.
*/
[[nodiscard]] std::vector<SRNdata> sr_matrix_elements(
  const std::vector<DiracSpinor> &orbs, DiracOperator::TensorOperator *h,
  MBPT::StructureRad *sr, ExternalField::CorePolarisation *dV,
  const SRNoptions &options, std::ostream &outstream = std::cout);

} // namespace Amplitudes
