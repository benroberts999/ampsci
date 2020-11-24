#pragma once
#include "CorePolarisation.hpp"
#include "TDHF.hpp"
#include <string>
#include <vector>
class Wavefunction;
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}

namespace HF {
class HartreeFock;
class Breit;
} // namespace HF

namespace ExternalField {

//! @brief
//! Similar to the time-dependent Hartree-Fock method, but expands dPsi
//! corrections using a basis to include core-polarisation (RPA) corrections to
//! matrix elements of some external field operator.

/*! @details
Like the TDHF method, but uses a basis:
Note: Benefit is we can include Breit into dV while using a basis (unlike
diagramRPA method)
*/
class TDHFbasis final : public TDHF {
public:
  TDHFbasis(const DiracOperator::TensorOperator *const h,
            const HF::HartreeFock *const hf,
            const std::vector<DiracSpinor> &basis);

private:
  // dPhi = X exp(-iwt) + Y exp(+iwt)
  // (H - e - w)X = -(h + dV - de)Phi
  // (H - e + w)Y = -(h* + dV* - de)Phi
  // X_c = sum_x X_x,
  // j(x)=j(c)-k,...,j(c)+k.  And: pi(x) = pi(c)*pi(h)
  // std::vector<std::vector<DiracSpinor>> m_X{};
  // std::vector<std::vector<DiracSpinor>> m_Y{};
  // can just write these to disk! Read them in, continue as per normal

  const std::vector<DiracSpinor> m_basis; // store copy?

public:
  //! @brief Solves TDHF equations self-consistantly for core electrons at
  //! frequency omega.
  //! @details Solves TDHF equations self-consistantly for core electrons at
  //! frequency omega. Will iterate up to a maximum of max_its. Set max_its=1
  //! to get first-order correction [note: no dampling is used for first
  //! itteration]. If print=true, will write progress to screen
  virtual void solve_core(const double omega, int max_its = 100,
                          const bool print = true) override final;

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
  //! st takes values: StateType::ket or StateType::bra
  DiracSpinor form_dPsi(const DiracSpinor &Fv, const double omega,
                        dPsiType XorY, const int kappa_beta,
                        const std::vector<DiracSpinor> &spectrum,
                        StateType st = StateType::ket) const;

  //! Forms \delta Psi_v for valence state Fv for all kappas (see solve_dPsi)
  std::vector<DiracSpinor> form_dPsis(const DiracSpinor &Fv, const double omega,
                                      dPsiType XorY,
                                      const std::vector<DiracSpinor> &spectrum,
                                      StateType st = StateType::ket) const;

public:
  TDHFbasis &operator=(const TDHFbasis &) = delete;
  TDHFbasis(const TDHFbasis &) = default;
  ~TDHFbasis() = default;
};

} // namespace ExternalField
