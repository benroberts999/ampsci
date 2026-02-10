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

enum class Method { TDHF, basis, diagram, none, Error };

inline Method ParseMethod(std::string_view str) {
  return qip::ci_compare(str, "TDHF")       ? Method::TDHF :
         qip::ci_compare(str, "true")       ? Method::TDHF :
         qip::ci_compare(str, "basis")      ? Method::basis :
         qip::ci_compare(str, "tdhf_basis") ? Method::basis :
         qip::ci_compare(str, "tdhfbasis")  ? Method::basis :
         qip::ci_compare(str, "diagram")    ? Method::diagram :
         qip::ci_compare(str, "rpad")       ? Method::diagram :
         qip::ci_compare(str, "rpa(d)")     ? Method::diagram :
         qip::ci_compare(str, "none")       ? Method::none :
         qip::ci_compare(str, "false")      ? Method::none :
         qip::ci_compare(str, "")           ? Method::none :
                                              Method::Error;
}

enum class dPsiType { X, Y };
enum class StateType { bra, ket }; // lhs, rhs

//! Virtual Core Polarisation class, for <a||dV||b>. See TDHF, DiagramRPA, etc.
class CorePolarisation {

protected:
  CorePolarisation(const DiracOperator::TensorOperator *const h)
      : m_h(h), m_rank(h->rank()), m_pi(h->parity()), m_imag(h->imaginaryQ()) {}

  CorePolarisation(int K, int Parity, bool imaginary)
      : m_h(nullptr), m_rank(K), m_pi(Parity), m_imag(imaginary) {}

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
  int rank() const { return m_rank; }
  int parity() const { return m_pi; }
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

  const DiracOperator::TensorOperator *get_operator() const { return m_h; }

  //! Returns RPA method
  virtual Method method() const = 0;

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
