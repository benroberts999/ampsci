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
  @brief Like @ref TDHF, but solves the TDHF equations via basis expansion.
  @details
  Same physics as @ref TDHF; see that class for full description.
  The perturbed core orbitals are found via the basis expansion
  \f[
    \varphi^a_\pm = \sum_n
      \frac{\ket{n}\matel{n}{t_\pm + \delta V_\pm}{a}}{\en_a - \en_n \pm \omega},
  \f]
  rather than by solving the inhomogeneous ODE directly.

  @warning Does not currently work correctly for frequency-dependent operators
  unless they depend only on the magnitude \f$ |\omega| \f$.
  The method assumes \f$ t_- = t_+^\dag \f$, whereas the correct relation is
  \f$ t_-(\omega) = t_+^\dag(-\omega) \f$.
  This will be fixed in a future update.
*/
class TDHFbasis final : public TDHF {
public:
  /*!
    @brief Constructs TDHFbasis for operator h using the provided basis.
    @param h     External field operator.
    @param hf    @ref HF::HartreeFock object defining the core.
    @param basis Single-particle basis used to expand \f$ \varphi^a_\pm \f$.
                 The entire basis is used; caller must ensure all required
                 states are present.
  */
  TDHFbasis(const DiracOperator::TensorOperator *const h,
            const HF::HartreeFock *const hf,
            const std::vector<DiracSpinor> &basis);

private:
  std::vector<DiracSpinor> m_basis{}; // store copy?

public:
  //! See @ref TDHF::solve_core(); same notes and warnings apply.
  virtual void solve_core(const double omega, int max_its = 100,
                          const bool print = true) override final;

  virtual Method method() const override final { return Method::basis; }

  /*!
    @brief Forms varphi^v_pm for valence state Fv: single kappa channel.
    @details
    Solves the TDHF equation via basis expansion for a single kappa channel;
    see @ref TDHF::solve_dPsi() for the equation. Returns \f$ \chi_\beta \f$
    for given kappa_beta, where
    \f[ X_{j,m} = (-1)^{j_\beta-m}tjs(j,k,j;-m,0,m)\chi_j. \f]

    @param Fv         Valence state \f$ \phi_v \f$.
    @param omega      Perturbation frequency \f$ \omega \f$.
    @param XorY       Selects X or Y solution; see @ref dPsiType.
    @param kappa_beta Kappa quantum number of the target channel.
    @param spectrum   Single-particle spectrum used to expand the solution.
    @param st         Bra or ket convention; see @ref StateType.
    @param incl_dV    Include the induced potential \f$ \delta V \f$ if true.

    @note To include correlations, use a basis with correlations.
    @note To exclude excitations to occupied states, remove them from the basis.
  */
  DiracSpinor form_dPsi(const DiracSpinor &Fv, const double omega,
                        dPsiType XorY, const int kappa_beta,
                        const std::vector<DiracSpinor> &spectrum,
                        StateType st = StateType::ket,
                        bool incl_dV = true) const;

  /*!
    @brief Forms varphi^v_pm for valence state Fv: all kappa channels.
    @details
    Calls @ref form_dPsi() for all allowed \f$ \kappa \f$ channels.

    @param Fv       Valence state \f$ \phi_v \f$.
    @param omega    Perturbation frequency \f$ \omega \f$.
    @param XorY     Selects the X or Y solution; see @ref dPsiType.
    @param spectrum Single-particle spectrum used to construct the solution.
    @param st       Bra or ket convention; see @ref StateType.
    @param incl_dV  Include the induced potential \f$ \delta V \f$ if true.
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
