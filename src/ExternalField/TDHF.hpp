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

/*!
  @brief Uses TDHF to include core-polarisation (RPA) corrections to matrix
  elements of an external field operator.

  @details
  Solves the set of TDHF equations
  \f[
    (H - \epsilon \pm \omega)\delta\psi_b = -(\delta V + \delta\epsilon_c)\psi_b
  \f]
  self-consistently for each electron in the core to determine \f$\delta V\f$.
  See 'ampsci.pdf' for a detailed physics description.

  \par Construction
  Requires a pointer to an operator (h) and a const @ref HF::HartreeFock object.

  \par Usage
  @ref solve_core (omega) solves the TDHF equations for a given frequency.
  Frequency should be positive (negative is allowed for testing only, with
  care). Can be re-run with a different frequency without restarting from
  scratch. @ref dV (Fa, Fb) then returns the correction to the matrix element:
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
  /*!
    @brief Constructs TDHF for operator h.

    @param h  External field operator.
    @param hf @ref HF::HartreeFock object defining the core.
  */
  TDHF(const DiracOperator::TensorOperator *const h,
       const HF::HartreeFock *const hf);

  /*!
    @brief Solves TDHF equations self-consistently for core electrons at
    frequency omega.

    @details Will iterate up to a maximum of max_its. Set max_its=1 to get
    first-order correction [note: no damping is used for first iteration].
    If print=true, will write progress to screen.
  */
  virtual void solve_core(double omega, int max_its = 100,
                          bool print = true) override;

  virtual Method method() const override { return Method::TDHF; }

  virtual void clear() override final;

  //! Returns reduced matrix element \f$\langle a||\delta V||b\rangle\f$,
  //! or the conjugate \f$\langle a||\delta V^\dagger||b\rangle\f$ if conj=true.
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb, bool conj) const;

  virtual double dV(const DiracSpinor &Fa,
                    const DiracSpinor &Fb) const override final;

  DiracSpinor dV_rhs(int kappa_n, const DiracSpinor &Fm,
                     bool conj = false) const override;

  //! Returns const ref to dPsi orbitals for given core orbital Fc.
  const std::vector<DiracSpinor> &get_dPsis(const DiracSpinor &Fc,
                                            dPsiType XorY) const;
  //! Returns const ref to dPsi orbital of given kappa.
  const DiracSpinor &get_dPsi_x(const DiracSpinor &Fc, dPsiType XorY,
                                const int kappa_x) const;

  /*! @brief
    Forms \f$\delta\psi_v\f$ for valence state Fv (including core pol.):
    single kappa.

    @details
    Solves
    \f[ (H + \Sigma - \epsilon - \omega)X
      = -(h + \delta V - \delta\epsilon)\psi \f]
    or
    \f[ (H + \Sigma - \epsilon + \omega)Y
      = -(h^\dagger + \delta V^\dagger - \delta\epsilon)\psi \f]
    Returns \f$ \chi_\beta \f$ for given kappa_beta, where
    \f[ X_{j,m} = (-1)^{j_\beta-m}tjs(j,k,j;-m,0,m)\chi_j \f]

    @param Fv         Valence state \f$\psi_v\f$.
    @param omega      Perturbation frequency \f$\omega\f$.
    @param XorY       Selects X or Y solution; see @ref dPsiType.
    @param kappa_beta Kappa quantum number of the target channel.
    @param Sigma      Optional correlation potential; see @ref MBPT::CorrelationPotential.
    @param st         Bra or ket convention; see @ref StateType.
    @param incl_dV    Include the induced potential \f$\delta V\f$ if true.
  */
  DiracSpinor
  solve_dPsi(const DiracSpinor &Fv, const double omega, dPsiType XorY,
             const int kappa_beta,
             const MBPT::CorrelationPotential *const Sigma = nullptr,
             StateType st = StateType::ket, bool incl_dV = true) const;

  //! Forms \f$\delta\psi_v\f$ for all kappa channels; see @ref solve_dPsi.
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
