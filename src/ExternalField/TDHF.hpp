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
  Solves the TDHF equations for each core orbital \f$ \phi_a \f$,
  \f[
    (h_{\rm HF} - \en_a \mp \omega)\varphi^a_\pm
      = -(t_\pm + \delta V_\pm - \delta\en^a_\pm)\phi_a,
  \f]
  self-consistently to determine \f$ \delta V_\pm \f$.
  See the @ref ExternalField namespace documentation for full physics description.

  Each core orbital acquires both a forward (\f$ \varphi^a_+ \f$, stored as X)
  and backward (\f$ \varphi^a_- \f$, stored as Y) correction. Both contribute
  to \f$ \delta V_\pm \f$:
  \f[
    \delta V_\pm \phi_i =
    \sum_a^{\rm core} \left[
      \matel{\phi_a}{Q}{\varphi^a_+}\phi_i - \matel{\phi_a}{Q}{\phi_i}\varphi^a_+
    + \matel{\varphi^a_-}{Q}{\phi_a}\phi_i - \matel{\varphi^a_-}{Q}{\phi_i}\phi_a
    \right].
  \f]
  It is via \f$ \delta V_\pm \f$ that the \f$ e^{\pm i\omega t} \f$ terms are mixed.

  \par Construction
  Requires a pointer to the forward operator (\f$ t_+ \f$), a const
  @ref HF::HartreeFock object, and an optional pointer to the backward
  operator (\f$ t_- \f$). If @p h_minus is omitted, \f$ t_- \f$ is taken
  to be the same pointer as \f$ t_+ \f$.
  The user is responsible for updating operator frequencies externally
  before calling solve_core().
  Note: Only need \f$t_-\f$ for operators that depend on the sign of the frequency (e.g., E1v). Most depend only on magnitude, in which case \f$t_- = t_+^\dag \f$, which is dealt with automatically.

  \par Usage
  @ref solve_core (omega) solves the TDHF equations for a given frequency.
  @ref dV (Fa, Fb) then returns the RPA correction to the reduced matrix element
  \f$ \redmatel{a}{\delta V}{b} \f$.
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
  std::vector<std::vector<DiracSpinor>> m_hFcore_minus{};
  // can just write these to disk! Read them in, continue as per normal

  const HF::HartreeFock *const p_hf;
  const std::vector<DiracSpinor> m_core;
  // const std::vector<double> m_Hmag;
  const double m_alpha;
  const HF::Breit *const p_VBr;
  // nb: m_h_plus := m_h is the one in CorePolarisation
  const DiracOperator::TensorOperator *const m_h_minus;

public:
  /*!
    @brief Constructs TDHF for operator h.

    @param h_plus  Forward operator \f$ t_+ \f$; must be set to positive
                   frequency before each call to solve_core().
    @param hf      @ref HF::HartreeFock object defining the core.
    @param h_minus Backward operator \f$ t_- \f$; if nullptr (default), uses
                   @p h_plus. Only needed when the operator is frequency-dependent
                   *and* depends on the **sign** of \f$ \omega \f$ (e.g. E1
                   velocity form). For operators that depend only on
                   \f$ |\omega| \f$ (e.g. M1), the default suffices.
  */
  TDHF(const DiracOperator::TensorOperator *const h_plus,
       const HF::HartreeFock *const hf,
       const DiracOperator::TensorOperator *const h_minus = nullptr);

  /*!
    @brief Solves TDHF equations self-consistently for core electrons at frequency omega.
    @details
    Iterates the TDHF equations until \f$ \delta V_\pm \f$ converges to within
    eps_target(), or @p max_its iterations have been performed.
    Can be re-run at a different frequency without restarting from scratch.

    @param omega    External-field frequency \f$ \omega \f$ in atomic units.
    @param max_its  Maximum number of iterations.
                   Set to 1 to get the first-order correction (no damping on
                   first iteration).
    @param print    If true, write convergence progress to screen.

    @note Frequency should be positive; negative is allowed but use with care.
    Unlike @ref DiagramRPA, TDHF solves for both \f$ \varphi^a_+ \f$ and
    \f$ \varphi^a_- \f$ simultaneously (see class docs), so a single call
    covers both contributions to \f$ \delta V_\pm \f$.

    ---

    @note Does not update the frequency of the operator itself; for frequency-dependent
    operators, update the operator frequency externally before calling.
  */
  virtual void solve_core(double omega, int max_its = 100,
                          bool print = true) override;

  virtual Method method() const override { return Method::TDHF; }

  virtual void clear() override final;

  //! Returns reduced matrix element \f$\redmatel{a}{\delta V}{b}\f$,
  //! or the conjugate \f$\redmatel{a}{\delta V^\dagger}{b}\f$ if conj=true.
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

  /*!
    @brief Forms \f$\varphi^v_\pm\f$ for valence state Fv (including core pol.):
    single kappa channel.
    @details
    Solves
    \f[ (h_{\rm HF} + \Sigma - \en_v - \omega)\varphi^v_+
      = -(t + \delta V - \delta\en^v)\phi_v \f]
    or
    \f[ (h_{\rm HF} + \Sigma - \en_v + \omega)\varphi^v_-
      = -(t^\dagger + \delta V^\dagger - \delta\en^v)\phi_v \f]
    Returns \f$ \chi_\beta \f$ for given kappa_beta, where
    \f[ X_{j,m} = (-1)^{j_\beta-m}tjs(j,k,j;-m,0,m)\chi_j \f]

    @param Fv         Valence state \f$\phi_v\f$.
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

  //! Forms \f$\varphi^v_\pm\f$ for all kappa channels; see @ref solve_dPsi.
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
  std::vector<std::vector<DiracSpinor>>
  form_hFcore(const DiracOperator::TensorOperator *h) const;
  // Solves the MS equations for all projections, single core state
  void solve_ms_core(std::vector<DiracSpinor> &dFb, const DiracSpinor &Fb,
                     const std::vector<DiracSpinor> &hFbs, const double omega,
                     dPsiType XorY, double eps_ms = 1.0e-9) const;

  /*!
    @brief Projects source `rhs` orthogonal to the near-resonant core states of
    a TDHF channel; returns the projected-out states and amplitudes.
    @details
    For the channel of angular momentum @p kappa_beta and energy
    \f$ e_0 = \en_b + ww \f$, the radial operator \f$ (h_l - e_0) \f$ is
    (near-)singular for components along any bound core state \f$ \phi_m \f$ of
    the same kappa with \f$ \en_m \approx e_0 \f$ -- exactly singular for the
    diagonal (\f$ \kappa_\beta = \kappa_b \f$, \f$ ww = 0 \f$) case, and
    near-singular for e.g. fine-structure partners. Those components cannot be
    found reliably from the Green's-function solve.

    Such states are identified by a relative nearness criterion
    (\f$ |e_0 - \en_m| < \eta\,|e_0 + \en_m| \f$, \f$ \eta \approx 0.1 \f$),
    removed from @p rhs, and returned so the caller can:
    
    (i) force the solution orthogonal to them during the inhomogeneous solve
    (keeping the remaining, non-resonant part well-conditioned), and
    
    (ii) restore the off-diagonal ones analytically via first-order PT,
    \f$ \matel{m}{dF}{} = \matel{m}{\rm src}{}/(\en_b \pm \omega - \en_m) \f$.
    
    The restored value equals the component the Green's-function solve would
    have produced naturally, so the result matches solving without projection
    but stays stable as the denominator becomes small. Only resonant states
    are touched: well-separated same-kappa core components are left in @p rhs
    and recovered naturally by the solve, since transitions into those occupied
    states cancel pairwise in \f$ \delta V \f$.

    @param rhs        Source spinor; modified in place (resonant part removed).
    @param core       Core orbitals.
    @param Fb         Core state \f$ \phi_b \f$ being perturbed.
    @param ww         Signed frequency (\f$ +\omega \f$ for X, \f$ -\omega \f$ for Y).
    @param kappa_beta Kappa of the solution channel.
    @return {states, amplitudes}: parallel vectors of the projected-out core
            states and the corresponding amplitudes \f$ \matel{m}{\rm src}{} \f$.
  */
  static std::pair<std::vector<const DiracSpinor *>, std::vector<double>>
  project_resonant(DiracSpinor &rhs, const std::vector<DiracSpinor> &core,
                   const DiracSpinor &Fb, double ww, int kappa_beta);

public:
  TDHF &operator=(const TDHF &) = delete;
  TDHF(const TDHF &) = default;
  ~TDHF() = default;
};

} // namespace ExternalField
