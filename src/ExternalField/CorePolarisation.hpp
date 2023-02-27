#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "qip/String.hpp"
#include <string>
#include <vector>
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}

//! Calculates many-body corrections (RPA) to matrix elements of external field
namespace ExternalField {

enum class method { TDHF, basis, diagram, none, Error };

inline method ParseMethod(std::string_view str) {
  return qip::ci_compare(str, "TDHF")    ? method::TDHF :
         qip::ci_compare(str, "basis")   ? method::basis :
         qip::ci_compare(str, "diagram") ? method::diagram :
         qip::ci_compare(str, "none")    ? method::none :
         qip::ci_compare(str, "")        ? method::none :
                                           method::Error;
}

enum class dPsiType { X, Y };
enum class StateType { bra, ket }; // lhs, rhs

//! Virtual Core Polarisation class, for <a||dV||b>. See TDHF, DiagramRPA
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

public:
  //! Returns eps (convergance) of last solve_core run
  double get_eps() const { return m_core_eps; }
  //! Returns omega (frequency) of last solve_core run
  double get_omega() const { return m_core_omega; }
  int get_rank() const { return m_rank; }
  int get_parity() const { return m_pi; }
  bool get_imagQ() const { return m_imag; }

  //! Solve RPA equations (for whichever method) for core.
  virtual void solve_core(const double omega, int max_its = 100,
                          const bool print = true) = 0;

  //! @brief Clears the dPsi orbitals (sets to zero)
  virtual void clear() = 0;

  //! @brief Calculate reduced matrix element <n||dV||m>
  virtual double dV(const DiracSpinor &Fn, const DiracSpinor &Fm) const = 0;

public:
  CorePolarisation &operator=(const CorePolarisation &) = delete;
  CorePolarisation(const CorePolarisation &) = default;
  virtual ~CorePolarisation() = default;
};

} // namespace ExternalField
