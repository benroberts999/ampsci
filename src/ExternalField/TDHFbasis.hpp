#pragma once
#include "CorePolarisation.hpp"
#include "TDHF.hpp"
#include <vector>
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}

namespace ExternalField {

/*!
  @brief Similar to TDHF, but uses a basis expansion to include
  core-polarisation (RPA) corrections to matrix elements of an external field.

  @details
  Like the @ref TDHF method, but uses a basis to expand the dPsi corrections.
  Solves
  \f[ (H + \delta V - \epsilon \pm \omega)\,\delta\psi = -(h + \delta V)\psi \f]
  by expanding \f$\delta\psi\f$ over the basis:
  \f[ \delta\psi = \sum_n \frac{|n\rangle\langle n|(h + \delta V)|\psi\rangle}{\epsilon - \epsilon_n \pm \omega} \f]
*/
class TDHFbasis final : public TDHF {
public:
  /*!
    @brief Constructs TDHFbasis for operator h using the provided basis.

    @param h     External field operator.
    @param hf    @ref HF::HartreeFock object defining the core.
    @param basis Single-particle basis used to expand dPsi. The entire basis
                 is used; it is the caller's responsibility to ensure all
                 required states are present.
  */
  TDHFbasis(const DiracOperator::TensorOperator *const h,
            const HF::HartreeFock *const hf,
            const std::vector<DiracSpinor> &basis);

private:
  std::vector<DiracSpinor> m_basis{}; // store copy?

public:
  virtual void solve_core(const double omega, int max_its = 100,
                          const bool print = true) override final;

  virtual Method method() const override final { return Method::basis; }

  /*! @brief
    Forms dF_v for valence state Fv (including core pol.): single kappa.

    @details
    Solves
    \f[ (H + \Sigma - \epsilon - \omega)X
      = -(h + \delta V - \delta\epsilon)\psi
    \f]
    or
    \f[
      (H + \Sigma - \epsilon + \omega)Y
        = -(h^\dagger + \delta V^\dagger - \delta\epsilon)\Psi
    \f]

    Returns \f$ \chi_\beta \f$ for given kappa_beta, where
    \f[
      X_{j,m} = (-1)^{j_\beta-m}tjs(j,k,j;-m,0,m)\chi_j
    \f]

    @param Fv         Valence state \f$\psi_v\f$.
    @param omega      Perturbation frequency \f$\omega\f$.
    @param XorY       Selects X or Y solution; see @ref dPsiType.
    @param kappa_beta Kappa quantum number of the target channel.
    @param spectrum   Single-particle spectrum used to expand the solution.
    @param st         Bra or ket convention; see @ref StateType.
    @param incl_dV    Include the induced potential \f$\delta V\f$ if true.

    @note To include correlations, use a basis with correlations.
    @note To exclude excitations to occupied states, remove them from the basis.
  */
  DiracSpinor form_dPsi(const DiracSpinor &Fv, const double omega,
                        dPsiType XorY, const int kappa_beta,
                        const std::vector<DiracSpinor> &spectrum,
                        StateType st = StateType::ket,
                        bool incl_dV = true) const;

  /*! @brief
    Forms the perturbed wavefunctions dF_v for a valence state for all allowed
    kappa channels.

    @details
    Solves the inhomogeneous equation for the perturbed wavefunction
    \f$\delta\psi_v\f$ of the valence state \f$\psi_v\f$ at frequency
    \f$\omega\f$. The solutions \f$[\delta\psi_v]_\kappa\f$ are returned for
    all allowed \f$\kappa\f$ values (see @ref form_dPsi).

    The calculation uses the supplied single-particle spectrum and may include
    the induced potential \f$\delta V\f$ if `incl_dV` is true.

    @param Fv       Valence state \f$\psi_v\f$.
    @param omega    Perturbation frequency \f$\omega\f$.
    @param XorY     Selects the X or Y solution; see @ref dPsiType.
    @param spectrum Single-particle spectrum used to construct the solution.
    @param st       Bra or ket convention; see @ref StateType.
    @param incl_dV  Include the induced potential \f$\delta V\f$ if true.
  */
  std::vector<DiracSpinor> form_dPsis(const DiracSpinor &Fv, const double omega,
                                      dPsiType XorY,
                                      const std::vector<DiracSpinor> &spectrum,
                                      StateType st = StateType::ket,
                                      bool incl_dV = true) const;

public:
  TDHFbasis &operator=(const TDHFbasis &) = delete;
  TDHFbasis(const TDHFbasis &) = default;
  ~TDHFbasis() = default;
};

} // namespace ExternalField
