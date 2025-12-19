#pragma once
#include "CorePolarisation.hpp"
#include "TDHF.hpp"
#include <vector>
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}

namespace ExternalField {

//! @brief
//! Similar to the time-dependent Hartree-Fock method, but expands dPsi
//! corrections using a basis to include core-polarisation (RPA) corrections to
//! matrix elements of some external field operator.

/*! @details
Like the TDHF method, but uses a basis:
Note: Benefit is we can include Breit into dV while using a basis (unlike
diagramRPA method). Downside, it's quite slow (though maybe could be made much
more efficient)
*/
class TDHFbasis final : public TDHF {
public:
  TDHFbasis(const DiracOperator::TensorOperator *const h,
            const HF::HartreeFock *const hf,
            const std::vector<DiracSpinor> &basis);

private:
  std::vector<DiracSpinor> m_basis{}; // store copy?

public:
  //! @brief Solves TDHF equations self-consistantly for core electrons at
  //! frequency omega, using basis expansion.
  //! @details Solves TDHF equations (using a basis to expand the states)
  //! self-consistantly for core electrons at frequency omega. Will iterate up
  //! to a maximum of max_its. Set max_its=1 to get first-order correction
  //! [note: no dampling is used for first itteration]. If print=true, will
  //! write progress to screen
  virtual void solve_core(const double omega, int max_its = 100,
                          const bool print = true) override final;

  //! Returns RPA method
  virtual Method method() const override final { return Method::basis; }

  //! Forms \delta Psi_v for valence state Fv (including core pol.) - 1 kappa
  //! @details
  //!  Solves
  //! \f[ (H + \Sigma - \epsilon - \omega)X = -(h + \delta V
  //! - \delta\epsilon)\psi \f]
  //! or
  //! \f[ (H + \Sigma - \epsilon + \omega)Y = -(h^\dagger + \delta V^\dagger
  //!   - \delta\epsilon)Psi\f]
  //! Returns \f$ \chi_\beta \f$ for given kappa_beta, where
  //! \f[ X_{j,m} = (-1)^{j_\beta-m}tjs(j,k,j;-m,0,m)\chi_j \f]
  //! XorY takes values: dPsiType::X or dPsiType::Y.
  //! st takes values: StateType::ket or StateType::bra.
  //! Solves by expanding over basis. To include correlations, use basis with
  //! correlations.
  DiracSpinor form_dPsi(const DiracSpinor &Fv, const double omega,
                        dPsiType XorY, const int kappa_beta,
                        const std::vector<DiracSpinor> &spectrum,
                        StateType st = StateType::ket,
                        bool incl_dV = true) const;

  //! Forms \delta Psi_v for valence state Fv for all kappas (see form_dPsi)
  std::vector<DiracSpinor> form_dPsis(const DiracSpinor &Fv, const double omega,
                                      dPsiType XorY,
                                      const std::vector<DiracSpinor> &spectrum,
                                      StateType st = StateType::ket,
                                      bool incl_dV = true) const;

public:
  TDHFbasis &operator=(const TDHFbasis &) = default;
  TDHFbasis(const TDHFbasis &) = default;
  ~TDHFbasis() = default;
};

} // namespace ExternalField
