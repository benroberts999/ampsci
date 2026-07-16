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
  //! channel data (homogeneous solutions, K amplitudes).
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
    @brief Anderson/Pulay (DIIS) acceleration of the outer self-consistency
    (default on).
    @details
    The TDHF fixed-point map is LINEAR in the corrections (dV is linear in
    {X,Y} and each channel solve is a linear inversion), so Anderson mixing
    (equivalent to preconditioned GMRES) converges wherever (1 - A) is
    non-singular -- including near the (auto)ionising resonances and the
    occupied-occupied near-degeneracies (omega ~ en_b - en_a), where the
    Picard multiplier exceeds 1 and DAMPED iteration diverges for any
    damping factor. Set false for the plain damped iteration (with staged
    Y; the bound-TDHF-style driver).
    @note Memory: the history holds 2 * depth (= 8) flattened copies of all
    the (X, Y) corrections, ~ 8 * n_channels * 2 * num_points doubles --
    of order a GB for a heavy atom on a dense grid, PER TDHFcntm instance.
    Callers parallelising over omega should account for this (the photoRPA
    module runs energies serially for this reason, with threads used inside
    each solve), or use set_anderson(false).
  */
  void set_anderson(bool anderson) { m_anderson = anderson; }

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
    @brief On-shell rescattering (unitarisation) of the standing-wave RRPA
    amplitudes: Kbar matrix and physical outgoing amplitudes.
    @details
    The driven RRPA solve uses real standing-wave boundary conditions
    (principal-value continuum), so the amplitudes D are the "mathematical"
    reaction-matrix ones: they have real poles wherever an RRPA eigenphase
    passes through pi/2 (e.g. just above a deep ionisation threshold, or
    through the Xe 4d giant resonance region), where the PHYSICAL
    (outgoing-wave) amplitude stays finite. Following the appendix of
    Johnson and Cheng, Phys. Rev. A 20, 978 (1979), we build the on-shell
    rescattering matrix from the N_P homogeneous coupled-channel solutions:
    seed each open channel j in turn with its regular continuum orbital and
    iterate the RRPA equations with NO external field; the converged
    channel-i amplitudes give
    \f[ \bar K_{ij} = \pi \matel{y^0_i}{R^{(j)}}{\,}, \f]
    (their Eq. A14). The physical outgoing amplitudes are then
    \f[ A = (1 - i\bar K)^{-1}\,\pi D , \f]
    which is finite through the standing-wave poles (D and Kbar share them),
    and per-channel \f$ |A_i|/\pi \f$ replaces \f$ |D_i| \f$ in the
    cross-section. Everything is computed in the phase reference of the
    local (KS-conditioned) continuum pair Freg/Firr -- the same reference
    the channel solves extract K against -- in which the driven solution
    has exactly zero regular component; |A_i| is reference-independent.

    @note Kbar is symmetric in exact arithmetic (reciprocity survives the
    diagonal phase rotation from the HF to the KS reference); `asymmetry`
    reports the worst relative deviation as a numerical diagnostic.
    @warning Requires a converged solve_core() at this omega. Each of the
    N_P seeded solves costs about as much as the driven solve.
  */
  struct Rescattering {
    //! Open channels, in matrix order: core-orbital index, channel kappa,
    //! photoelectron energy en = en_a + omega
    struct Channel {
      std::size_t i_core;
      int kappa;
      double en;
    };
    std::vector<Channel> channels{};
    //! On-shell rescattering matrix Kbar_ij (KS phase reference)
    LinAlg::Matrix<double> Kbar{};
    //! Driven standing-wave amplitudes D_i = K_i/pi (KS phase reference)
    std::vector<double> D{};
    //! Physical amplitudes |A_i|/pi, A = (1 - i*Kbar)^{-1} pi*D
    std::vector<double> D_phys{};
    //! Worst |Kbar_ij - Kbar_ji| relative to the largest |Kbar| element
    double asymmetry{0.0};
  };
  /*!
    @brief See @ref Rescattering. max_its/eps as the driven solve_core.
    @param max_its   Maximum iterations for each seeded solve.
    @param print     Print each seeded solve's convergence (labelled by the
                     seeded channel). Forces the seeds to run serially.
    @param parallel  Solve the seed columns in parallel (one seeded solve
                     per thread; each thread holds its own copy of the
                     corrections plus an Anderson history -- similar memory
                     per thread to an omega-parallel driven solve). Ignored
                     when @p print is set.
  */
  Rescattering rescattering(int max_its = 40, bool print = false,
                            bool parallel = false) const;

  /*!
    @brief Unitarised (physical) amplitude |A|/pi for the open channel of
    hole Fa with photoelectron kappa_e; replaces |D| in the cross-section.
    @details
    The physical amplitudes are computed for ALL open channels in one go
    (the rescattering solves couple every channel; see @ref rescattering);
    the result is cached, so the first call after solve_core() does the
    work and later calls are lookups. Returns 0 for closed (or excluded)
    channels.
    @note This is a magnitude: the physical amplitude is complex, so the
    sign of the real standing-wave amplitude has no meaning here. Fine for
    cross-sections (squares).
  */
  double D_phys(const DiracSpinor &Fa, int kappa_e);

  //! Kbar asymmetry diagnostic of the cached rescattering (see @ref
  //! Rescattering); 0.0 if not yet computed (no D_phys call).
  double rescattering_asymmetry() const {
    return m_resc ? m_resc->asymmetry : 0.0;
  }

  //! Worst |K - pi*D| / K_max over the open channels of the last
  //! solve_core() (K_max = largest channel amplitude): the K = pi*D
  //! identity holds to ~1% for healthy channels, so a gross violation
  //! flags an unreliable channel solve (marginal grid resolution) even
  //! when the SCF itself converged. Scaled to the DOMINANT amplitude, so
  //! weak/zero-crossing channels cannot fire the alarm on a negligible
  //! absolute error. A warning is printed (if print) when it exceeds ~5%.
  double KpiD_dev() const { return m_KpiD; }

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
      across iterations, omegas, and the seeded (rescattering) solves;
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
  // by build_ycc() on the first solve_core(), shared by copies (seeded
  // rescattering solves run on copies).
  using YccTable = std::vector<std::vector<std::vector<double>>>;
  std::shared_ptr<const YccTable> m_ycc{};
  void build_ycc();
  const std::vector<double> &ycc_get(std::size_t i, std::size_t j, int l) const;
  // Continuum option flags (see setters above)
  bool m_suppress_open{false};
  bool m_staged_Y{true};
  bool m_anderson{true};

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
  // Cached rescattering result for D_phys()/rescattering_asymmetry();
  // invalidated by solve_core() and clear(). max_its of the last solve_core
  // is reused for the seeded solves.
  std::optional<Rescattering> m_resc{};
  int m_max_its{40};
  // Worst K = pi*D violation of the last solve_core (see KpiD_dev())
  double m_KpiD{0.0};

  // Homogeneous (seeded, no external field) solve state -- set only inside
  // rescattering(), on a copy of the driven-solved object. The seeded
  // channel's m_X entry stores the TOTAL (seed + correction), so dV sees the
  // seed; its own solve is for the correction, with the seed's exchange
  // deficit (vnl - U_KS)*seed added to the source. See rescattering().
  struct SeedInfo {
    std::size_t ib;      // core-orbital index of the seeded channel
    std::size_t be;      // channel index (within m_X[ib]) of the seed
    DiracSpinor seed;    // the seed: cached Freg of that channel
    DiracSpinor deficit; // (vexFa - vexFa_1el - U_KS) applied to the seed
  };
  std::optional<SeedInfo> m_seed{};

  // Runs the seeded homogeneous TDHF (no external field) for open channel
  // (ib, be): sets m_seed, zeroes X/Y, inserts the seed, and iterates with
  // the Anderson driver. On convergence the per-channel K amplitudes in
  // m_ch hold column (ib, be) of the rescattering matrix. Mutates *this
  // (call on a copy).
  void solve_homogeneous(std::size_t ib, std::size_t be, int max_its,
                         bool print);

  // (Re)builds m_ch for this omega: openness flags, and the homogeneous
  // continuum pair Freg/Firr for each open channel (barely-open channels
  // handled by the outward-extension pair construction; high-energy
  // channels truncated at Freg.max_pt()). If print, warns on the backstop
  // exclusion (solveContinuum returned zero: grid resolves nothing).
  void prepare_channels(double omega, bool print);

  // Drivers for solve_core (after the shared set-up): plain damped (staged)
  // fixed-point iteration, and Anderson/DIIS-accelerated iteration (see
  // set_anderson). warm_start: continuing at the same omega (first
  // iteration must be damped).
  void solve_core_damped(double omega, int max_its, bool print,
                         bool warm_start);
  void solve_core_anderson(double omega, int max_its, bool print);

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
  // hFb may be nullptr: no external-field source (homogeneous/seeded mode);
  // if ch is the seeded channel (see m_seed), the solve is for the
  // correction (seed deficit added to the source; total stored back).
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
