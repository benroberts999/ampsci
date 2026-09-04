#pragma once
#include "Coulomb/YkTable.hpp"
#include "TDHF.hpp"
#include <complex>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ExternalField {

/*!
  @brief Continuum TDHF/RRPA with outgoing-wave (complex) boundary
  conditions: core polarisation for any external field above ionisation
  threshold(s), giving the physical amplitudes from ONE driven solve.

  @details
  As @ref TDHF, but supports open (continuum) X/+ channels, where a core
  orbital is ionised: en_+ = en_a + omega > 0. Those channels are solved
  with the OUTGOING-wave boundary condition, so the corrections are
  complex, and the converged outgoing amplitude of each open channel is
  directly the physical (rescattering-included) amplitude. This is the
  alternative to the standing-wave solve plus on-shell K matrix of
  @ref TDHFcntm (Johnson method): one complex SCF (two real channel solves
  per channel per iteration) replaces N_open + 1 real SCFs, and no K
  matrix is ever built. The two routes give the same amplitudes (same
  convention, same channel order, same phase reference) and cross-validate
  each other.

  \par Physics scheme (shared with TDHFcntm)
  - Ionised orbitals are treated in the V^{N-1} (residual ion) potential
    via an exact rearrangement: the one-electron self-interaction
    V^a_0 = y^0_aa + X_a of the hole is moved from the lagged dV source
    into the static Hamiltonian of each OPEN channel and compensated in the
    source (@ref hole_compensation), so the fixed point is the unchanged
    TDHF one but the channel is solved with the correct Z_ion = 1 Coulomb
    tail. Bound channels (closed orbitals, and every Y/- partner) stay in
    plain V^N form, as bound TDHF, solved with Anderson mixing
    (@ref solveMixedState_cntm).
  - Occupied components are kept in the open channels (their dV
    contributions cancel pairwise across channels); only the diagonal
    phi_a is projected, for even-parity operators (kappa_beta = kappa_a),
    together with the norm-conservation term de*phi_a in the source, so
    the response to t -> identity is exactly zero.

  \par Method
  In the local-KS channel reference of the cached pair (Freg ~ cos, Firr ~
  sin, so the outgoing wave is Freg + i Firr), the standing-wave solve
  returns phi_P ~ K Firr with K = pi <Freg|S>; the outgoing solution of the
  same equation is phi_+ = phi_P - i K Freg ~ -iK (Freg + i Firr), i.e.
  G^+ = G^P - i pi |Freg><Freg|. With the nonlocal exchange iterated in the
  source, let F^HF = Freg + (standing response to its exchange deficit) be
  the exchange-dressed regular solution of the channel, F^HF ~ Freg + K_ex
  Firr (cached per open channel per omega). For a complex source
  S = S_r + i S_i the outgoing exchange-dressed solution is
  \f[ \varphi_+ = \varphi_P[S_r] + i\varphi_P[S_i] - i K_+ F^{HF},
      \qquad K_+ = \frac{K_r + iK_i}{1 + iK_{ex}} , \f]
  two ordinary standing-wave solves (@ref solveContinuumMixedState) plus a
  multiple of the cached F^HF.

  The TDHF map itself needs no change: with the Y/- sets storing the
  complex CONJUGATE of the true Y corrections, the real-linear dV builder
  applied separately to the real and imaginary sets is the exact complex
  TDHF source (delta V_+ is linear in X and in conj(Y); the external field
  is real in the code's convention and drives the real part only). The
  real and imaginary sets couple solely through the outgoing addition in
  the open channels.

  At the fixed point K_+ = (1 + i Kbar)^{-1} pi D, the complex conjugate of
  the Johnson amplitude A = (1 - i Kbar)^{-1} pi D (outgoing-wave response
  vs incoming-wave-normalised final state; |A| and every same-channel
  relative phase are identical). @ref A_phys returns conj(K_+)/pi.

  \par Usage
  As TDHF: @ref solve_core (omega), then either
  - @ref dV_complex (Fe, Fa), with Fe the energy-normalised V^{N-1} HF
    continuum state of hole Fa (ContinuumOrbitals::solveContinuumHF with
    hole_particle and force_orthog): the complex dressed correction, so
    that |<Fe|t|a> + dV_complex(Fe, Fa)| is the physical amplitude (phase
    in the HF reference; accuracy limited by Fe's inner region); or
  - @ref A_phys / @ref D_phys (Fa, kappa_e): the physical amplitude read
    off the converged solution directly, with no continuum bra at all
    (local-KS phase reference).
  Below all thresholds every channel is bound, the imaginary sets stay
  zero, and dV matches bound TDHF. One instance serves many operators
  (@ref set_operator; the omega caches are kept within a (rank, parity)
  class); @ref prepare does the omega-only work on a master instance whose
  copies serve per-thread driven solves.
*/
class TDHFcomplex : public TDHF {

public:
  /*!
    @brief Constructs the outgoing-wave continuum TDHF for operator h; see
    @ref TDHF::TDHF.
    @details The core-only tables (the core-core Coulomb screening functions
    of the batched dV builder) are built here, once; copies share them.
  */
  TDHFcomplex(const DiracOperator::TensorOperator *const h_plus,
              const HF::HartreeFock *const hf,
              const DiracOperator::TensorOperator *const h_minus = nullptr);

  //! Open-channel label: core-orbital index, channel kappa, and the
  //! ionised-electron energy en = en_a + omega
  struct Channel {
    std::size_t i_core;
    int kappa;
    double en;
  };

  /*!
    @brief Solves the complex (outgoing-wave) continuum TDHF equations
    self-consistently at frequency omega; see class description.
    @param omega    External-field frequency (atomic units).
    @param max_its  Maximum number of (Anderson) iterations; 1 gives the
                    first-order (outgoing) correction.
    @param print    If true, write convergence progress to screen.
  */
  void solve_core(double omega, int max_its = 100, bool print = true) override;

  //! Clears the corrections (as TDHF::clear) and the cached continuum
  //! channel data
  void clear() override;

  /*!
    @brief Switches the instance to a different external-field operator
    (with its conjugate partner, as in the constructor).
    @details The core-only tables are kept always. The omega-only caches
    (channel openness, continuum pairs, exchange-dressed incident waves)
    depend on the operator only through its rank and parity, so they are
    kept when those are unchanged (e.g. the two E1 gauges; every multipole
    of one (rank, parity) class). The previous corrections are discarded.
    Updating an operator's frequency or momentum transfer IN PLACE needs
    no call here: the next solve_core() warm-starts from the previous
    corrections.
  */
  void set_operator(const DiracOperator::TensorOperator *h_plus,
                    const DiracOperator::TensorOperator *h_minus = nullptr);

  /*!
    @brief The operator-independent work at omega: channel openness, the
    continuum pair of each open channel, and its exchange-dressed incident
    wave F^HF with K_ex. No driven solve.
    @details Cheap (one channel solve per open channel). solve_core() calls
    it; call it directly on a master instance whose copies (carrying the
    caches) do the driven solves of many operators / momentum transfers.
    Nothing is rebuilt while omega, rank, and parity are unchanged.
  */
  void prepare(double omega);

  //! Open (usable) channels at the prepared omega, in the order of
  //! @ref A_phys. Empty before prepare() or solve_core().
  std::vector<Channel> channel_list() const;

  //! If true, the open (continuum) X/+ channels are excluded (held at zero)
  //! rather than solved. Diagnostic option; default false.
  void set_suppress_open(bool suppress_open) {
    m_suppress_open = suppress_open;
  }

  /*!
    @brief Physical amplitudes A/pi of every open channel of the converged
    solve, as complex numbers, in channel_list() order.
    @details A = conj(K_+) (see class description): the Johnson convention
    A = (1 - i Kbar)^{-1} pi D, in the local-KS (conditioning-potential)
    phase reference of each channel. The physical outgoing-wave amplitude
    carries a further channel phase exp(i(delta_i + sigma_i)) (HF and
    Coulomb phase shifts) that is NOT included, so the RELATIVE phase of
    two operators' amplitudes in the SAME channel is physical (e.g. the
    interference term Re(A_t conj(A_L))), while the absolute phase, and
    the phase between different channels, are not.
    @warning Requires a converged solve_core() at this omega.
  */
  std::vector<std::complex<double>> A_phys() const;

  //! Physical amplitude A/pi (see @ref A_phys) of the open channel of hole
  //! Fa with ionised-electron kappa_e; 0 for closed (or excluded) channels.
  std::complex<double> A_phys(const DiracSpinor &Fa, int kappa_e) const;

  //! |A_phys| for the channel of hole Fa with ionised-electron kappa_e: the
  //! magnitude that enters cross-sections.
  double D_phys(const DiracSpinor &Fa, int kappa_e) const;

  /*!
    @brief Per-open-channel results for hole orbital Fa after solve_core():
    the complex outgoing amplitude K_+ and the internal amplitude
    D = <Freg|S>, S the total (complex) effective source of the converged
    solve (physical source, hole compensation, iterated exchange
    remainder), with the SAME local-KS Freg the solver used.
    @details The radial solve satisfies K_+ = pi*D identically, so the
    residual is a stringent check of the solve and the boundary-condition
    extraction (see @ref KpiD_dev).
  */
  struct OutgoingChannel {
    int kappa;
    double en;
    std::complex<double> K;
    std::complex<double> D;
  };
  //! See @ref OutgoingChannel. Empty if solve_core() not yet run, or no open
  //! channels for Fa at this omega.
  std::vector<OutgoingChannel> outgoing_channels(const DiracSpinor &Fa) const;

  /*!
    @brief Reduced ME of the (complex) induced potential, <a||dV||b>, or the
    conjugate <a||dV^dagger||b> if en_b > en_a; as @ref TDHF::dV, complex.
    @details
    Bound a and b: real (and equal to bound TDHF) below every threshold,
    complex in general above.

    Continuum bra (en_a > 0): @p Fa must be the energy-normalised V^{N-1}
    continuum state of hole @p Fb (ContinuumOrbitals::solveContinuumHF with
    hole_particle and force_orthog), and the source carries the hole
    compensation V^b_0 phi (see @ref hole_compensation), matching the
    rearranged channel equations: <Fe || dV + V^b_0 phi || Fb>. Then
    <Fe|t|b> + dV_complex(Fe, Fb) is the outgoing-wave (physical) dressed
    amplitude: its modulus is |A| of @ref A_phys and its phase is in the
    HF reference (K_+ exp(i delta_ex)). Accuracy is limited by Fe's inner
    region; @ref A_phys needs no bra.
    @warning Requires solve_core() at the matching omega, and (continuum
    bra) Fe.kappa() to be one of the channels of Fb.
  */
  std::complex<double> dV_complex(const DiracSpinor &Fa,
                                  const DiracSpinor &Fb) const;

  //! Real part of @ref dV_complex, for bound a and b (exact below every
  //! threshold: the class then acts as bound TDHF). Asserts on a continuum
  //! state: use dV_complex.
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb) const override;
  using TDHF::dV;

  //! Worst |K_+ - pi*D| / K_max over the open channels of the last
  //! solve_core() (K_max = largest channel amplitude): the identity holds
  //! to ~1% for healthy channels, so a gross violation flags an unreliable
  //! channel solve (marginal grid resolution) even when the SCF converged.
  //! Scaled to the DOMINANT amplitude, so weak channels cannot fire the
  //! alarm on a negligible absolute error. Never printed: query it.
  double KpiD_dev() const { return m_KpiD; }

  //! Label ("shell,kappa") of the channel attaining KpiD_dev(); empty if
  //! there are no open channels.
  const std::string &KpiD_worst_channel() const { return m_KpiD_lab; }

  //! Open channels at the prepared omega dropped from the response because
  //! the continuum solve returned zero (the grid resolves nothing at that
  //! energy); they are zeroed. Empty when none were dropped (normal).
  const std::string &excluded_channels() const { return m_excluded; }

  //! Sum of the squared norms of the imaginary parts of all corrections
  //! (X and Y); exactly zero when no channel is open. Diagnostic.
  double imaginary_norm2() const;

  /*!
    @brief The hole (self-interaction) compensation term V^a_0 chi =
    (y^0_aa + X_a) chi for ionised orbital Fa.
    @details Compensates, in the source, the one-electron self-interaction
    moved into the static V^{N-1} Hamiltonian of the open channels (vl -
    y^0_aa on the local side, the one-electron self-exchange via the
    solver's Fhole option). Its +y^0_aa*chi piece cancels pointwise the 1/r
    tail hidden in the b=a part of dV*phi_a, leaving a short-ranged total
    source.
  */
  DiracSpinor hole_compensation(const DiracSpinor &Fa,
                                const DiracSpinor &chi) const;

  //! Batched dV sources for one TDHF iteration; see @ref dV_rhs_all.
  struct dVrhs_XY {
    //! X[ib][be] = dV_rhs(kappa_be, core[ib], conj = false); indexed as m_X
    std::vector<std::vector<DiracSpinor>> X{};
    //! Y[ib][be] = dV_rhs(kappa_be, core[ib], conj = true); empty if skipped
    std::vector<std::vector<DiracSpinor>> Y{};
  };
  /*!
    @brief Builds [dV phi_a]_beta for EVERY (core orbital x channel x X/Y)
    task of one TDHF iteration from the correction sets X, Y (indexed as
    m_X, m_Y) in a single batched pass; same result as TDHF::dV_rhs per
    task.
    @details The per-task dV_rhs recomputes every radial Coulomb (yk)
    screening function from scratch. Batched, they are shared: the direct
    (Q) parts of both W terms collapse into ONE screening function per task
    type, built once per call; the exchange (P) part of the first W term
    needs only the fixed core-core y^l(Fb, Fa) (built once per instance);
    the exchange part of the second W term needs y^l(eta_beta, Fa), once
    per core orbital. Called twice per iteration here: on the real sets
    and on the imaginary sets (the builder is real-linear).
    @param include_Y  Also build the conjugate (Y) sources.
    @note Falls back to per-task TDHF::dV_rhs when a Breit potential is
    present (the Breit dV term is not restructured). Public for testing.
  */
  dVrhs_XY dV_rhs_all(const std::vector<std::vector<DiracSpinor>> &X,
                      const std::vector<std::vector<DiracSpinor>> &Y,
                      bool include_Y) const;

private:
  // Core-core Coulomb screening functions y^l(core_i, core_j), used by the
  // exchange part of dV_rhs_all(). Fixed for the life of the instance, so
  // built on construction; shared_ptr so copies share it.
  std::shared_ptr<const Coulomb::YkTable> m_ycc;

  // Imaginary parts of the corrections, indexed as m_X / m_Y (which hold
  // the real parts). The Y sets store conj(Y) (see class description).
  std::vector<std::vector<DiracSpinor>> m_Xi{};
  std::vector<std::vector<DiracSpinor>> m_Yi{};

  // Per-(core orbital x channel) continuum data, cached at fixed omega:
  // openness; the homogeneous pair Freg/Firr at en_+ (built once, in the
  // fixed conditioning potential vl - y^0_aa + U_KS); the exchange-dressed
  // incident wave F^HF ~ Freg + K_ex Firr; and the complex outgoing
  // amplitude K_+ of the latest solve. Indexed as m_X.
  struct ContinuumChannel {
    bool open{false};
    DiracSpinor Freg;
    DiracSpinor Firr;
    DiracSpinor Fhf;
    double K_ex{0.0};
    std::complex<double> K{0.0};
  };
  std::vector<std::vector<ContinuumChannel>> m_channels{};
  // The omega the caches were built at (rebuild when it changes)
  double m_omega{-1.0};
  // A driven solution at m_omega is held
  bool m_have_solution{false};
  bool m_suppress_open{false};
  // Worst K = pi*D violation of the last solve_core (see KpiD_dev()), and
  // the channel that attained it
  double m_KpiD{0.0};
  std::string m_KpiD_lab{};
  // Channels dropped by the zero-pair backstop at m_omega
  std::string m_excluded{};

  // (Re)shapes the imaginary sets to m_X / m_Y, zeroed
  void zero_imaginary_sets();

  // (Re)builds m_channels for this omega: openness flags; for each open
  // channel the pair Freg/Firr (barely-open channels by the outward-
  // extension construction; high-energy channels truncated at
  // Freg.max_pt()), then F^HF and K_ex. Channels whose continuum solve
  // returned zero are recorded in m_excluded and dispatched as zero.
  void prepare_channels(double omega);

  // Open (usable) channels in fixed (core orbital, channel) order; fills
  // *channels and returns the (i_core, i_channel) index list
  std::vector<std::pair<std::size_t, std::size_t>>
  list_open_channels(std::vector<Channel> *channels) const;

  // (i_core, i_channel) of the channel of core orbital Fb with kappa
  std::pair<std::size_t, std::size_t> channel_index(const DiracSpinor &Fb,
                                                    int kappa) const;

  // Anderson/Pulay (DIIS) accelerated self-consistency driver: the TDHF
  // fixed-point map is linear in the corrections, so DIIS (preconditioned
  // GMRES on the doubled, real + imaginary, state) converges wherever
  // (1 - A) is non-singular, including through resonance windows. State =
  // all (X, Y, Xi, Yi) corrections and the complex channel amplitudes.
  void solve_core_anderson(double omega, int max_its, bool print);

  // Single undamped TDHF iteration, all (core orbital x channel x X/Y)
  // solves task-flattened; returns the convergence measure
  std::pair<double, std::string> tdhf_core_it_complex(double omega);

  // Outgoing-wave solve of one open X/+ channel: real and imaginary parts
  // of the previous iterate in / new iterate out; K_+ written to the
  // channel cache. Sources, diagonal projection, and lagged hole
  // compensation for both parts, then the outgoing addition.
  void solve_channel_outgoing(DiracSpinor *X_re, DiracSpinor *X_im,
                              ContinuumChannel *channel, const DiracSpinor &Fb,
                              const DiracSpinor &hFb, const DiracSpinor &dV_re,
                              const DiracSpinor &dV_im, double omega,
                              double eps_ms) const;

  // Bound channel (closed orbital, or the Y/- partner of an ionised
  // orbital): plain V^N bound solve (Anderson), as bound TDHF. hFb may be
  // nullptr (no external-field source: the imaginary part)
  void solve_channel_bound(DiracSpinor *dF_beta, const DiracSpinor &Fb,
                           const DiracSpinor *hFb, const DiracSpinor &dV_src,
                           double omega, dPsiType XorY, double eps_ms) const;

  // Convergence measure: bound channels the relative L2 change of
  // (X_re, X_im) over the norm of ALL X channels (physically-negligible
  // bound remnants must not gate the SCF); open channels
  // |dK_+|^2 / max(|K_+|^2, floor) (box-independent, physical). K_old
  // indexed as m_channels.
  std::pair<double, std::string> eps_complex(
    const std::vector<std::vector<DiracSpinor>> &Xs_re,
    const std::vector<std::vector<DiracSpinor>> &Xs_im,
    const std::vector<std::vector<std::complex<double>>> &K_old) const;

public:
  TDHFcomplex &operator=(const TDHFcomplex &) = delete;
  TDHFcomplex(const TDHFcomplex &) = default;
  ~TDHFcomplex() = default;
};

} // namespace ExternalField
