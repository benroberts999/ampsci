#pragma once
#include "LinAlg/Matrix.hpp"
#include "TDHF.hpp"
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace ExternalField {

/*!
  @brief Continuum (photoionisation) TDHF/RRPA: core polarisation at
  frequencies above ionisation threshold(s), following the Johnson method.

  @details
  As @ref TDHF, but supports open (continuum) X/+ channels, where a core
  orbital is ionised: en_+ = en_a + omega > 0. All continuum logic lives in
  this class; the bound TDHF base class is untouched.

  Two changes relative to the bound solver:

  1. Ionised orbitals (en_a + omega > 0) are treated in the V^{N-1}
  (residual ion) potential, via an exact rearrangement of the TDHF
  equations: the ONE-electron (spherically averaged) self-interaction
  V^a_0 = y^0_aa + X_a of the hole orbital (X_a the one-electron
  self-exchange, sign included as in vexFa: X_a chi = vexFa_1el(chi))
  is moved from the lagged dV
  source into the static Hamiltonian,
  \f[ (h_{\rm HF} - V^a_0 - \en_\pm)\varphi_\pm
      = -(t + \delta V - \delta\en)\phi_a - V^a_0\bar\varphi_\pm , \f]
  so the fixed point is the unchanged TDHF one, but the X/+ channel can be
  solved with the standing-wave continuum solve
  (@ref solveContinuumMixedState) in a potential with the correct
  residual-ion (Z_ion = 1) Coulomb tail. The rearrangement is applied per
  channel, to the OPEN channels only: every bound channel (closed orbitals,
  and the bound Y/- partners of ionised orbitals) stays in plain V^N form,
  exactly as bound TDHF. The compensating (lagged) source
  term is built in ONE place, @ref hole_compensation, so the hole-particle
  treatment can be adjusted/verified in isolation.
  Note: only ONE electron is removed, not the subshell: the dynamic
  response of the other [j_a]-1 same-subshell electrons remains in dV
  (cf. Johnson et al., Phys. Scr. 21, 409 (1980), Sec 3.1).

  2. Bound channels (the Y/- partners and all channels of closed orbitals)
  are solved with Anderson/Pulay mixing (@ref solveMixedState_cntm) rather
  than damped iteration: at photoionisation frequencies a bound channel's
  en_0 = en_a + omega can sit near a spurious eigenvalue of the local
  preconditioner, where damped iteration diverges.

  Even-parity operators (diagonal channels, kappa_beta = kappa_a; e.g. the
  K=0 temporal multipole t^0 = j0(qr) of electron-impact ionisation): the
  norm-conservation
  (Lagrange-multiplier) term -de*phi_a, de = <a|(t + dV)phi_a>, is
  subtracted from the physical source of each open diagonal solve (the V^a_0
  compensation is part of the operator rearrangement and is NOT projected),
  and the solution is kept orthogonal to phi_a. This is the exact
  counterpart of the bound solver's conditioning projection; in particular
  the response to t -> identity (j0(qr) at q -> 0) is exactly zero.

  \par Usage
  As TDHF: @ref solve_core (omega) at a frequency (typically above at least
  one ionisation threshold; below all thresholds every channel is bound and
  the result matches bound TDHF). Matrix elements with a continuum final
  state are then evaluated with @ref dV_cntm (not @ref dV), which applies
  the matching V^{N-1} source term.
*/
class TDHFcntm : public TDHF {

public:
  /*!
    @brief Constructs continuum TDHF for operator h; see @ref TDHF::TDHF.
  */
  TDHFcntm(const DiracOperator::TensorOperator *const h_plus,
           const HF::HartreeFock *const hf,
           const DiracOperator::TensorOperator *const h_minus = nullptr);

  /*!
    @brief Solves the (continuum) TDHF equations self-consistently at
    frequency omega; see class description and @ref TDHF::solve_core.
    @param omega    External-field frequency (atomic units); typically above
                    (at least one) ionisation threshold.
    @param max_its  Maximum number of iterations. Set to 1 for the
                    first-order correction.
    @param print    If true, write convergence progress to screen.
  */
  virtual void solve_core(double omega, int max_its = 100,
                          bool print = true) override;

  //! Clears the dPsi corrections (as TDHF::clear) and the cached continuum
  //! channel data (continuum pairs, K amplitudes).
  virtual void clear() override;

  //! If true, the open (continuum) X/+ channels are excluded (held at zero)
  //! rather than solved: dV then contains only the closed (bound) part of
  //! the core response. Diagnostic option; default false.
  void set_suppress_open(bool suppress_open) {
    m_suppress_open = suppress_open;
  }

  //! Staged iteration (default on; damped mode only): after the first-order
  //! step, the Y/- corrections are held frozen while X is iterated, and
  //! released once X has roughly converged. Improves stability, per Johnson
  //! et al. (1979 appendix): "it is convenient to ignore the
  //! negative-frequency orbitals entirely at the first stages".
  void set_staged_Y(bool staged_Y) { m_staged_Y = staged_Y; }

  /*!
    @brief Per-open-channel results for hole orbital Fa, after solve_core().
    @details
    For each open channel (kappa, en_+ = en_a + omega > 0): the standing-wave
    (K-matrix) amplitude K of the converged correction, and the internal
    amplitude D = <F_reg|S>, with the SAME (local-KS) F_reg the solver used
    and S the total effective source of the converged solve (including the
    iterated nonlocal-exchange remainder). The radial solve satisfies
    K = pi*D identically, so the residual is a stringent check of the
    solve, the boundary-condition extraction, and the short-rangedness of
    the source (a long-ranged residual makes both sides box-dependent).
    @note The production amplitude uses the full HF V^{N-1} continuum bra
    (ContinuumOrbitals) and dV_cntm(); the D here is a diagnostic only.
  */
  struct OpenChannel {
    int kappa;
    double en;
    double K;
    double D;
  };
  //! See @ref OpenChannel. Empty if solve_core() not yet run, or no open
  //! channels for Fa at this omega.
  std::vector<OpenChannel> open_channels(const DiracSpinor &Fa) const;

  /*!
    @brief Reduced ME of dV for a continuum final state, consistent with the
    V^{N-1} treatment of the photoelectron: <Fe || dV + V^a_0 phi || Fa>.
    @details
    As @ref TDHF::dV, but adds the one-electron (spherically averaged)
    self-interaction term V^a_0 phi_pm = (y^0_aa + X_a) phi_pm to the
    source, matching the rearranged continuum TDHF equations (see class
    description): the static V^a_0 lives in the (directly inverted) V^{N-1}
    Hamiltonian, compensated by this lagged source term, so the fixed point
    is the exact TDHF one. The +y^0_aa phi piece cancels (pointwise) the
    1/r tail hidden in the b=a part of dV*phi_a, making the
    continuum-continuum matrix element box-independent.

    Use this for photoionisation amplitudes, with @p Fe the energy-normalised
    V^{N-1} continuum state of hole @p Fa (including the exchange part, see
    @ref ContinuumOrbitals::solveContinuumHF with subtract_self); use
    @ref dV for bound-bound.

    @note The diagonal de (norm-conservation) term of the TDHF source is NOT
    included: it multiplies phi_a, and @p Fe is orthogonalised against the
    core, so it contributes nothing to the matrix element.

    @warning Requires solve_core() to have been run at the matching omega,
    and @p Fe.kappa() to be one of the dPsi channels of @p Fa (guaranteed
    when the operator's selection rules connect Fe and Fa).
  */
  double dV_cntm(const DiracSpinor &Fe, const DiracSpinor &Fa) const;

  /*!
    @brief Source term [dV' phi_a]_beta for an IONISED orbital Fa: dV_rhs
    plus the hole compensation. The single hole-particle seam.
    @details
    \f[ [\delta V'\phi_a]_\beta = [\delta V\phi_a]_\beta
        + V^a_0\,\chi_\beta , \f]
    where chi is the (lagged) correction spinor for this channel (previous
    iterate during the TDHF; converged channel for the dressed amplitude).
    Public so tests can probe the short-rangedness of the total source.
    @note Does NOT include the diagonal de (norm-conservation) projection of
    the TDHF source (even-parity operators): that applies to the physical
    (t + dV) part only, and is done inside the channel solve.
  */
  DiracSpinor dVprime_rhs(int kappa_beta, const DiracSpinor &Fa, bool conj,
                          const DiracSpinor &chi) const;

  /*!
    @brief The hole (self-interaction) compensation term V^a_0 chi =
    (y^0_aa + X_a) chi for ionised orbital Fa.
    @details
    Compensates, in the source, the one-electron self-interaction moved into
    the static V^{N-1} Hamiltonian (vl - y^0_aa on the local side, and the
    one-electron self-exchange via the solver's Fhole option). Its
    +y^0_aa*chi piece must cancel pointwise the 1/r tail hidden in the b=a
    part of dV*phi_a, leaving a short-ranged total source; this is verified
    numerically in the [cntmrpa] tests.
  */
  DiracSpinor hole_compensation(const DiracSpinor &Fa,
                                const DiracSpinor &chi) const;

  //! Batched dV sources for one TDHF iteration; see @ref dV_rhs_all.
  struct DVRhsAll {
    //! X[ib][be] = dV_rhs(kappa_be, core[ib], conj = false); indexed as m_X
    std::vector<std::vector<DiracSpinor>> X{};
    //! Y[ib][be] = dV_rhs(kappa_be, core[ib], conj = true); empty if skipped
    std::vector<std::vector<DiracSpinor>> Y{};
  };
  /*!
    @brief Builds [dV phi_a]_beta for EVERY (core orbital x channel x X/Y)
    task of one TDHF iteration in a single batched pass; same result as
    calling TDHF::dV_rhs per task.
    @details
    The per-task dV_rhs recomputes every radial Coulomb (yk) screening
    function from scratch, which dominates the TDHF iteration cost. Batched,
    the yk functions are shared:
    - the direct (Q) parts of both W terms collapse into ONE radial
      screening function per task type (X/Y), built once per iteration
      (the (-1)^((j_n - j_beta)) phase separates into per-channel factors);
    - the exchange (P) part of the first W term needs only the core-core
      y^l(Fb, Fa): fixed functions, computed once per instance and shared
      across iterations and omegas;
    - the exchange (P) part of the second W term needs y^l(eta_beta, Fa):
      built once per core orbital a, shared across its target channels.
    @param include_Y  Also build the conjugate (Y) sources; skipped while Y
                      is frozen (staged iteration).
    @note Falls back to per-task TDHF::dV_rhs when a Breit potential is
    present (the Breit dV term is not restructured). Public for testing.
  */
  DVRhsAll dV_rhs_all(bool include_Y) const;

private:
  // Core-core Coulomb screening functions y^l(core_i, core_j), for every
  // pair (j <= i, index i*(i+1)/2 + j) and every multipole l allowed by the
  // Ck selection rules (index (l - lmin)/2, lmin from Angular::kminmax_Ck).
  // Fixed for the life of the instance (the core never changes); built once
  // by build_ycc() on the first solve_core(), shared by copies.
  using YccTable = std::vector<std::vector<std::vector<double>>>;
  std::shared_ptr<const YccTable> m_ycc{};
  void build_ycc();
  const std::vector<double> &ycc_get(std::size_t i, std::size_t j, int l) const;
  // Continuum option flags (see setters above)
  bool m_suppress_open{false};
  bool m_staged_Y{true};

  // Per-(core orbital x channel) continuum data, cached at fixed omega:
  // openness, the homogeneous pair Freg/Firr at en_+ (built ONCE, in the
  // fixed conditioning potential vl - y^0_aa + U_KS), and the standing-wave
  // K amplitude from the latest solve. Indexed as m_X.
  struct CntmChannel {
    bool open{false};
    double K{0.0};
    DiracSpinor Freg;
    DiracSpinor Firr;
    CntmChannel(bool t_open, DiracSpinor t_Freg, DiracSpinor t_Firr)
      : open(t_open), Freg(std::move(t_Freg)), Firr(std::move(t_Firr)) {}
  };
  std::vector<std::vector<CntmChannel>> m_ch{};
  // The omega the caches were built at (rebuild when it changes)
  double m_omega{-1.0};
  // (Re)builds m_ch for this omega: openness flags (open AND resolvable on
  // the radial box), and the homogeneous continuum pair Freg/Firr for each
  // open channel. If print, warns when open channels are excluded as
  // unresolvable (en_+ too close to that shell's threshold for the box).
  void prepare_channels(double omega, bool print);

  // Driver for solve_core (after the shared set-up): plain damped (staged)
  // fixed-point iteration. warm_start: continuing at the same omega (first
  // iteration must be damped).
  void solve_core_damped(double omega, int max_its, bool print,
                         bool warm_start);

  // Single TDHF iteration, all (core orbital x channel x X/Y) solves
  // task-flattened; continuum-aware dispatch. Mirrors TDHF::tdhf_core_it.
  // If !include_Y, the Y/- solves are skipped (Y frozen; staged iteration).
  std::pair<double, std::string>
  tdhf_core_it_cntm(double omega, double eta_damp, bool include_Y);

  // Convergence measure: bound channels contribute the relative L2 change
  // of the (undamped) X spinors (as TDHF::eps_dPsi); open channels the
  // relative change of the standing-wave K amplitude (the L2 norm of a
  // box-truncated continuum spinor is box-dependent; K is physical).
  // K_old holds the previous iteration's K, indexed as m_ch.
  std::pair<double, std::string>
  eps_cntm(const std::vector<std::vector<DiracSpinor>> &Xnew,
           const std::vector<std::vector<double>> &K_old) const;

  // Single-channel solve with the continuum (V^{N-1}) dispatch: open X/+
  // channels via the standing-wave continuum (forward) solve (Freg/Firr
  // reused from *ch; K written back to it; ch may be nullptr for Y), the
  // bound Y/- partner of an ionised orbital with the matching hole-particle
  // term, closed orbitals with the plain bound (Anderson) solve.
  // dV_src: this task's [dV phi_a]_beta, from dV_rhs_all() (built once per
  // iteration; replaces the per-task dV_rhs call).
  void solve_channel_cntm(DiracSpinor *dF_beta, CntmChannel *ch,
                          const DiracSpinor &Fb, const DiracSpinor *hFb,
                          const DiracSpinor &dV_src, double omega,
                          dPsiType XorY, double eps_ms) const;

public:
  TDHFcntm &operator=(const TDHFcntm &) = delete;
  TDHFcntm(const TDHFcntm &) = default;
  ~TDHFcntm() = default;
};

} // namespace ExternalField
