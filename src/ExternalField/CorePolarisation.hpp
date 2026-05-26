#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "qip/String.hpp"
#include <cassert>
#include <string>
#include <vector>
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}

/*!
  @brief Core-polarisation (RPA) corrections to matrix elements of an external field.
  @details
  Provides classes and functions for computing all-order core-polarisation
  corrections to matrix elements of an external field operator.

  In the presence of an external field of frequency \f$ \omega \f$, the
  time-dependent operator is
  \f[
    T(t) = t_+ e^{-i\omega t} + t_- e^{+i\omega t},
  \f]
  where \f$ t_+ = t^k_q(\omega) \f$ is an irreducible tensor operator of rank \f$ k \f$
  and \f$ t_- = t_+^\dag(-\omega) \f$.
  Each orbital, including those in the core, acquires a first-order perturbation,
  \f[
    \delta\phi_a(t) = \varphi^a_+ e^{-i\omega t} + \varphi^a_- e^{+i\omega t}.
  \f]
  Since the core orbitals are perturbed, the Hartree-Fock potential is also perturbed:
  \f[
    \delta V_\pm \phi_i =
    \sum_a^{\rm core} \left[
      \matel{\phi_a}{Q}{\varphi^a_+}\phi_i - \matel{\phi_a}{Q}{\phi_i}\varphi^a_+
    + \matel{\varphi^a_-}{Q}{\phi_a}\phi_i - \matel{\varphi^a_-}{Q}{\phi_i}\phi_a
    \right].
  \f]

  The resulting 
  core-polarisation corrections to matrix elements are given by
  \f[
    \matel{b}{t_\pm}{a} \to \matel{b}{t_\pm + \delta V_\pm}{a},
  \f]
  where \f$ \delta V_\pm \f$ is the correction to the HF potential arising
  from the perturbed core orbitals \f$ \{\varphi^a_\pm\} \f$.

  Since \f$ \delta V \f$ is solved self-consistently, this accounts for core
  polarisation to all orders in the Coulomb interaction.

  Two equivalent methods are implemented:

  - TDHF (@ref TDHF, @ref TDHFbasis): solves the TDHF equations
  \f[
    (h_{\rm HF} - \en_a \mp \omega)\varphi^a_\pm
      = -(t_\pm + \delta V_\pm - \delta\en^a_\pm)\phi_a
  \f]
  self-consistently for all core orbitals. The corrections can also be
  found via the basis expansion
  \f[
    \varphi^a_\pm = \sum_n \frac{\ket{n}\matel{n}{t_\pm + \delta V_\pm}{a}}
                              {\en_a - \en_n \pm \omega}.
  \f]
  See Dzuba et al. (1984).

  - Diagram/Goldstone (@ref DiagramRPA): evaluates the four core-polarisation
  diagrams directly,
  \f[
    \redmatel{w}{\delta V_\pm}{v} =
    \sum_{na} \frac{(-1)^{k+w-n}}{[k]} \left(
      \frac{T_{na}\,W^k_{wavn}}{\en_a - \en_n \pm \omega}
      + \frac{(-1)^{a-n} T_{an}\,W^k_{wnva}}{\en_a - \en_n \mp \omega}
    \right),
  \f]
  with \f$ T_{ij} = t_{ij} + \redmatel{i}{\delta V_\pm}{j} \f$ iterated to
  convergence. See Johnson et al., Phys. Rev. A 21, 409 (1980).

  Use @ref make_rpa to construct the appropriate object from a method string,
  and @ref calcMatrixElements to compute matrix elements including the correction.
*/
namespace ExternalField {

//! Available RPA/core-polarisation methods
enum class Method { TDHF, basis, diagram, none, Error };

//! Parses method string to Method enum (case-insensitive)
inline Method ParseMethod(std::string_view str) {
  return qip::ci_compare(str, "TDHF")       ? Method::TDHF :
         qip::ci_compare(str, "true")       ? Method::TDHF :
         qip::ci_compare(str, "default")    ? Method::TDHF :
         qip::ci_compare(str, "basis")      ? Method::basis :
         qip::ci_compare(str, "tdhf_basis") ? Method::basis :
         qip::ci_compare(str, "tdhfbasis")  ? Method::basis :
         qip::ci_compare(str, "diagram")    ? Method::diagram :
         qip::ci_compare(str, "diagramRPA") ? Method::diagram :
         qip::ci_compare(str, "rpad")       ? Method::diagram :
         qip::ci_compare(str, "rpa(d)")     ? Method::diagram :
         qip::ci_compare(str, "none")       ? Method::none :
         qip::ci_compare(str, "false")      ? Method::none :
         qip::ci_compare(str, "")           ? Method::none :
                                              Method::Error;
}

/*!
  @brief Selects the perturbed orbital: X = varphi_+, Y = varphi_-.
  @details
  Corresponds to the two first-order corrections to a core orbital,
  \f[ \delta\phi_a(t) = \varphi^a_+ e^{-i\omega t} + \varphi^a_- e^{+i\omega t}. \f]
  X selects \f$ \varphi^a_+ \f$ (absorption/forward), Y selects
  \f$ \varphi^a_- \f$ (emission/backward).
*/
enum class dPsiType { X, Y };
//! Whether the state is a bra or ket
enum class StateType { bra, ket }; // lhs, rhs

/*!
  @brief Virtual base class for core-polarisation (RPA); computes dV corrections.
  @details
  Defines the interface for all RPA/core-polarisation methods. Concrete
  implementations are @ref TDHF, @ref TDHFbasis, and @ref DiagramRPA.
  See the @ref ExternalField namespace documentation for notation and physics.

  @note
  Stores a raw pointer to the external-field operator @p h passed at
  construction. That operator must remain alive for the lifetime of this object.

  ---

  @note
  For frequency-dependent operators, updating the operator frequency externally
  will affect results. @ref solve_core() should be re-called after any such update.

  ---

  @note
  Calling @ref solve_core() with a different freuqnecy will _only_ chnage the frequency used
  to solve the TDHF/RPA equations. It wil **not** change the frequency of the operator itself.
  For frequency-dependent operators, you must update the freuqency of the operator first.
  See @ref DiracOperator::TensorOperator .
  Since this class stores just a pointer to the operator, that's all you need to do.
*/
class CorePolarisation {

protected:
  CorePolarisation(const DiracOperator::TensorOperator *const h)
    : m_h(h), m_rank(h->rank()), m_pi(h->parity()), m_imag(h->imaginaryQ()) {}

protected:
  const DiracOperator::TensorOperator *m_h;
  double m_core_eps{1.0};
  int m_core_its{0};
  double m_core_omega{0.0};
  int m_rank;
  int m_pi;
  bool m_imag;

  double m_eta{0.4};
  double m_eps{1.0e-10};

public:
  //! Returns eps (convergance) of last solve_core run
  double last_eps() const { return m_core_eps; }
  //! Returns its (# of iterations) of last solve_core run
  double last_its() const { return m_core_its; }
  //! Returns omega (frequency) of last solve_core run
  double last_omega() const { return m_core_omega; }
  //! Rank of the operator
  int rank() const { return m_rank; }
  //! Parity of the operator
  int parity() const { return m_pi; }
  //! Returns true if the operator is imaginary
  bool imagQ() const { return m_imag; }

  //! Convergance target
  double &eps_target() { return m_eps; }
  //! Convergance target
  double eps_target() const { return m_eps; }

  //! Damping factor; 0 means no damping. Must have 0 <= eta < 1
  double eta() const { return m_eta; }
  //! Set/update damping factor; 0 means no damping. Must have 0 <= eta < 1
  void set_eta(double eta) {
    assert(eta >= 0.0 && eta < 1 && "Must have 0 <= eta < 1");
    m_eta = eta;
  }

  //! Returns RPA method
  virtual Method method() const = 0;

  /*!
    @brief Solve for delta_V_pm self-consistently for all core orbitals at frequency omega.
    @details
    Iterates the RPA equations (TDHF or diagram, depending on the implementation)
    until the correction \f$ \delta V_\pm(\omega) \f$ converges to within
    eps_target(), or @p max_its iterations are reached.

    @param omega    External-field frequency \f$ \omega \f$ in atomic units.
                   May be zero (static limit) or negative.
    @param max_its  Maximum number of iterations.
                   - 0: no iterations; dV() returns 0.
                   - 1: single iteration; dV() returns lowest-order correction.
    @param print    If true, print convergence progress to stdout.

    @note Does not update the frequency of the operator itself; for frequency-dependent
    operators, update the operator frequency externally before calling.
  */
  virtual void solve_core(double omega, int max_its = 100,
                          bool print = true) = 0;

  //! @brief Clears the internal state back to pre solve_core()
  virtual void clear() = 0;

  //! @brief Returns reduced matrix element <n||dV_pm||m> (see namespace doc for dV_pm)
  virtual double dV(const DiracSpinor &Fn, const DiracSpinor &Fm) const = 0;

  //! @brief Returns [dV_pm * phi_m]_kappa: RHS of TDHF eq., projected onto kappa (see namespace doc)
  virtual DiracSpinor dV_rhs(int kappa, const DiracSpinor &Fm,
                             bool conj = false) const {
    // XXX Remove this implementation (make pure virtual) once j_L killed
    (void)kappa;
    (void)Fm;
    (void)conj;
    assert(false && "This should be made pure virtual");
    return Fm;
  }

public:
  CorePolarisation &operator=(const CorePolarisation &) = delete;
  CorePolarisation(const CorePolarisation &) = default;
  virtual ~CorePolarisation() = default;
};

} // namespace ExternalField
