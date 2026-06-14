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

  // If true, eps is the relative change |dPsi| (sqrt of the ratio); if false,
  // the squared ratio. See eps_dPsi().
  bool m_eps_sqrt{false};

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

  /*!
    @brief Solves the TDHF equations at a frequency above an ionisation
    threshold, where core orbitals are ionised (en_a + omega > 0): the
    continuum RPA for photoionisation.
    @details
    As @ref solve_core, but supports open (continuum) X/+ channels.
    Two changes relative to solve_core:

    1. Ionised orbitals (en_a + omega > 0) are treated in the V^{N-1}
    (residual ion) potential, via an exact rearrangement of the TDHF
    equations: the ONE-electron (spherically averaged) self-interaction
    V^a_0 = y^0_aa - X_a of the hole orbital is moved from the lagged dV
    source into the static Hamiltonian,
    \f[ (h_{\rm HF} - V^a_0 - \en_\pm)\varphi_\pm
        = -(t + \delta V - \delta\en)\phi_a - V^a_0\bar\varphi_\pm , \f]
    so the fixed point is the unchanged TDHF one, but the X/+ channel can be
    solved with the standing-wave continuum solve
    (@ref solveContinuumMixedState) in a potential with the correct
    residual-ion (Z_ion = 1) Coulomb tail. The Y/- partner stays bound.
    Note: only ONE electron is removed, not the subshell: the dynamic
    response of the other [j_a]-1 same-subshell electrons remains in dV.

    2. The self-consistency iteration uses Pulay/DIIS extrapolation rather
    than damped iteration: the continuum fixed-point map has spectral radius
    of order 1, where damped iteration cannot converge.

    @param omega    External-field frequency (atomic units). Should be above
                    (at least one) ionisation threshold; below all thresholds
                    every orbital is treated as in solve_core.
    @param max_its  Maximum number of iterations. Set to 1 for the
                    first-order correction.
    @param print    If true, write convergence progress to screen.
    @param suppress_open  If true, the open (continuum) X/+ channels are
                    excluded (held at zero) rather than solved: dV then
                    contains only the closed (bound) part of the core
                    response, using only the standard bound machinery.

    @note Matrix elements with a continuum final state must be evaluated with
    @ref dV_cont (not @ref dV), which applies the matching V^{N-1} source term.
  */
  void solve_core_cntm(double omega, int max_its = 100, bool print = true,
                       bool suppress_open = false);

  virtual Method method() const override { return Method::TDHF; }

  virtual void clear() override final;

  //! Returns reduced matrix element \f$\redmatel{a}{\delta V}{b}\f$,
  //! or the conjugate \f$\redmatel{a}{\delta V^\dagger}{b}\f$ if conj=true.
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb, bool conj) const;

  virtual double dV(const DiracSpinor &Fa,
                    const DiracSpinor &Fb) const override final;

  DiracSpinor dV_rhs(int kappa_n, const DiracSpinor &Fm,
                     bool conj = false) const override;

  /*!
    @brief Reduced ME of dV for a continuum final state, consistent with the
    V^{N-1} treatment of the photoelectron: <Fe || dV + V^a_0 phi || Fa>.
    @details
    As @ref dV, but adds the one-electron (spherically averaged)
    self-interaction term V^a_0 phi_pm = (y^0_aa - X_a) phi_pm to the source,
    matching the rearranged continuum TDHF equations (see
    @ref solve_core_cntm): the static V^a_0 lives in the (directly inverted)
    V^{N-1} Hamiltonian, compensated by this lagged source term, so the fixed
    point is the exact TDHF one. The +y^0_aa phi piece cancels (pointwise) the
    1/r tail hidden in the b=a part of dV*phi_a, making the
    continuum-continuum matrix element box-independent.

    Use this for photoionisation amplitudes, with @p Fe the energy-normalised
    V^{N-1} continuum state of hole @p Fa (including the exchange part, see
    @ref ContinuumOrbitals::solveContinuumHF with subtract_self); use @ref dV
    for bound-bound.

    @warning Requires solve_core_cntm() to have been run at the matching
    omega, and @p Fe.kappa() to be one of the dPsi channels of @p Fa
    (guaranteed when the operator's selection rules connect Fe and Fa).
  */
  double dV_cont(const DiracSpinor &Fe, const DiracSpinor &Fa) const;

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
  // As solve_ms_core(), but a single channel (one projection). Thread-safe;
  // used to parallelise tdhf_core_it() over (orbital x channel x X/Y).
  void solve_ms_core_b(DiracSpinor &dF_beta, const DiracSpinor &Fb,
                       const DiracSpinor &hFb, const double omega,
                       dPsiType XorY, double eps_ms = 1.0e-9) const;

  // Continuum (photoionisation) version of solve_ms_core: handles ionised
  // orbitals (en_a + omega > 0) in the V^{N-1} potential -- open X/+ channels
  // via the standing-wave continuum (forward) solve, the bound Y/- partner
  // with the matching hole-particle term. Used only by solve_core_cntm()
  // (dispatched from solve_ms_core when m_cntm_mode is set).
  void solve_ms_core_cntm(std::vector<DiracSpinor> &dFb, const DiracSpinor &Fb,
                          const std::vector<DiracSpinor> &hFbs,
                          const double omega, dPsiType XorY,
                          double eps_ms = 1.0e-9) const;

  // Continuum (photoionisation) mode flags: set only by solve_core_cntm()
  // for its duration. Both false => solve_ms_core behaviour is identical to
  // the standard bound TDHF.
  bool m_cntm_mode{false};
  bool m_suppress_open{false};

  // Convergence (eps): the relative L2 change of the (undamped) X spinors,
  // Sum|dX|^2 / Sum|X_new|^2 (summed over all channels); returns its sqrt
  // -- the relative change -- if @p relative. Returns {eps, worst-channel}.
  std::pair<double, std::string>
  eps_dPsi(const std::vector<std::vector<DiracSpinor>> &Xnew,
           bool relative) const;

public:
  TDHF &operator=(const TDHF &) = delete;
  TDHF(const TDHF &) = default;
  ~TDHF() = default;
};

} // namespace ExternalField
