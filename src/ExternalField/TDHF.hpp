#pragma once
#include "CorePolarisation.hpp"
#include <string>
#include <vector>
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace MBPT {
class CorrelationPotential;
}

namespace HF {
class HartreeFock;
class Breit;
} // namespace HF

namespace ExternalField {

//! Uses time-dependent Hartree-Fock method to include core-polarisation
//! (RPA) corrections to matrix elements of some external field operator.
/*! @details
Solves set of TDHF equations
\f[ 
  (H -\epsilon \pm \omega)\delta\psi_b = -(\delta V  + \delta\epsilon_c)\psi_b
\f] 
self consistantly for each electron in the core to determine dV. (See
'ampsci.pdf' for detailed physics description). There is an option to limit
the maximum number of iterations; set to 1 to get the first-order correction
(nb: no damping is used for first iteration).
\par Construction
Requires a pointer to an operator (h), a const HF object (const pointer).
\par Usage
solve_core(omega) solves TDHF eqs for given frequency. Frequency should be
positive, but is allowed to be negative (use as a test only, with care). Can be
run again with a different frequency, typically does not need to be re-started
from scratch. Then, dV(Fa,Fb) returns the correction to the matrix element:
\f[ \langle \phi_a || \delta V || \phi_b \rangle \f]
*/
class TDHF : public CorePolarisation {

protected:
  // dPhi = X exp(-iwt) + Y exp(+iwt)
  // (H - e - w)X = -(h + dV - de)Phi
  // (H - e + w)Y = -(h* + dV* - de)Phi
  // X_c = sum_x X_x,
  // j(x)=j(c)-k,...,j(c)+k.  And: pi(x) = pi(c)*pi(h)
  std::vector<std::vector<DiracSpinor>> m_X{};
  std::vector<std::vector<DiracSpinor>> m_Y{};
  std::vector<std::vector<DiracSpinor>> m_hFcore{};
  // can just write these to disk! Read them in, continue as per normal

  const HF::HartreeFock *const p_hf;
  const std::vector<DiracSpinor> m_core;
  // const std::vector<double> m_Hmag;
  const double m_alpha;
  const HF::Breit *const p_VBr;

public:
  //! Contruct TDHF operator: takes pointer to operator and to HF object
  TDHF(const DiracOperator::TensorOperator *const h,
       const HF::HartreeFock *const hf);

  //! Solves TDHF equations self-consistantly for core electrons at frequency
  //! omega.
  /*! @details Will iterate up to a maximum of max_its. Set max_its=1
     to get first-order correction [note: no damping is used for first
     itteration]. If print=true, will write progress to screen.
  */
  virtual void solve_core(const double omega, int max_its = 100,
                          const bool print = true) override;

  //! Returns RPA method
  virtual Method method() const override { return Method::TDHF; }

  //! Clears the dPsi orbitals (sets to zero)
  virtual void clear() override final;

  //! Calculate reduced matrix element <a||dV||b> or <a||dV*||b>.
  //! Will exclude orbital 'Fexcl' from sum over core (for tests only)
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb, bool conj) const;

  //! As above, but automatically determines if 'conjugate' version
  //! required (Based on sign of [en_a-en_b])
  virtual double dV(const DiracSpinor &Fa,
                    const DiracSpinor &Fb) const override final;

  //! Returns "reduced partial matrix element RHS": dV||Fb}.
  //! Note: Fa * dV_rhs(..) equiv to dV(..)
  DiracSpinor dV_rhs(const int kappa_n, const DiracSpinor &Fm,
                     bool conj = false) const;

  //! Returns const ref to dPsi orbitals for given core orbital Fc
  const std::vector<DiracSpinor> &get_dPsis(const DiracSpinor &Fc,
                                            dPsiType XorY) const;
  //! Returns const reference to dPsi orbital of given kappa
  const DiracSpinor &get_dPsi_x(const DiracSpinor &Fc, dPsiType XorY,
                                const int kappa_x) const;

  //! Forms \delta Psi_v for valence state Fv (including core pol.) - 1 kappa
  /*! @details 
   Solves
   \f[ (H + \Sigma - \epsilon - \omega)X = -(h + \delta V
   - \delta\epsilon)\psi \f]
   or
   \f[ (H + \Sigma - \epsilon + \omega)Y = -(h^\dagger + \delta V^\dagger
     - \delta\epsilon)Psi\f]
   Returns \f$ \chi_\beta \f$ for given kappa_beta, where
   \f[ X_{j,m} = (-1)^{j_\beta-m}tjs(j,k,j;-m,0,m)\chi_j \f]
   XorY takes values: dPsiType::X or dPsiType::Y.
   st takes values: StateType::ket or StateType::bra
  */
  DiracSpinor
  solve_dPsi(const DiracSpinor &Fv, const double omega, dPsiType XorY,
             const int kappa_beta,
             const MBPT::CorrelationPotential *const Sigma = nullptr,
             StateType st = StateType::ket, bool incl_dV = true) const;
  //! Forms \delta Psi_v for valence state Fv for all kappas (see solve_dPsi)
  std::vector<DiracSpinor>
  solve_dPsis(const DiracSpinor &Fv, const double omega, dPsiType XorY,
              const MBPT::CorrelationPotential *const Sigma = nullptr,
              StateType st = StateType::ket, bool incl_dV = true) const;

  // //! Writes dPsi (f-component) to textfile
  // void print(const std::string &ofname = "dPsi.txt") const;

private:
  void initialise_dPsi();

  // Single iteration of TDHF equations
  std::pair<double, std::string> tdhf_core_it(double omega, double eta_damp);
  // Forms set of h*Fc for all core orbitals and all projections
  std::vector<std::vector<DiracSpinor>> form_hFcore() const;
  // Solves the MS equations for all projections, single core state
  void solve_ms_core(std::vector<DiracSpinor> &dFb, const DiracSpinor &Fb,
                     const std::vector<DiracSpinor> &hFbs, const double omega,
                     dPsiType XorY, double eps_ms = 1.0e-9) const;

public:
  TDHF &operator=(const TDHF &) = delete;
  TDHF(const TDHF &) = default;
  ~TDHF() = default;
};

} // namespace ExternalField
