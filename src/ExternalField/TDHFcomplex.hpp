#pragma once
#include "TDHF.hpp"
#include <complex>
#include <string>
#include <utility>
#include <vector>

namespace ExternalField {

/*!
  @brief TDHF (RPA) core polarisation above ionisation threshold(s): the open
  channels are solved with the outgoing-wave boundary condition, so the
  corrections and the corrected matrix elements are complex.

  @details
  As @ref TDHF, for a frequency above the ionisation threshold of one or more
  core orbitals, 
  
  \f[ 
    \en^a_+ = \en_a + \omega > 0. 
  \f] 
  
  Every X (+) channel of such an orbital is open, and describes the ejected 
  electron. For it, the solution of the TDHF equation regular at the origin is 
  not unique (the continuum orbital at en_+ may be added freely); 
  the physical solution is the purely outgoing wave at large r, 
  which is complex. 
  
  All other channels (closed orbitals, and every Y (-) partner) are bound, 
  and are solved as in TDHF.
  Through dV every correction is complex: the real and imaginary parts
  are stored as separate sets (m_X, m_Y and m_Xi, m_Yi) and iterated
  together. Below every threshold the imaginary sets are exactly zero and the
  class reproduces TDHF.

  In dV the (-) corrections enter as complex conjugates, so the Y sets store
  \f$ (\varphi_-)^* \f$ and the (real-linear) dV builder of TDHF applies to
  each set unchanged; the external field drives the real set only, and the
  two sets couple solely through the outgoing boundary condition.

  \par Open channels
  The ejected electron moves in the field of the residual ion: the
  self-interaction of one electron in the ionised orbital,
  \f$ V^a_0 = y^0_{aa} + X_a \f$ (direct and one-electron exchange), is moved
  from the HF Hamiltonian into the source. This is an exact rearrangement:
  the channel Hamiltonian has the correct -1/r tail, and the source is short
  ranged (the y^0_aa term cancels, pointwise, the 1/r tail of the a term of
  dV phi_a).

  With the real regular and irregular solutions F_reg, F_irr of the local
  channel Hamiltonian at en_+ (F_reg + i F_irr is outgoing), the standing-wave
  solve returns \f$ \varphi_P \to K F_{\rm irr} \f$, and the outgoing solution
  is \f$ \varphi_+ = \varphi_P - iK F_{\rm reg} \f$. 
  
  With the non-local
  exchange iterated in the source, F_reg is replaced by the exchange-dressed
  regular solution \f$ F^{\rm HF} \to F_{\rm reg} + K_{\rm ex} F_{\rm irr} \f$
  (built once per channel per omega), and
  
  \f[ K_+ = \frac{K}{1 + iK_{\rm ex}}, \f]
  
  and
  
  \f[
      \varphi_+ = \varphi_P - iK_+F^{\rm HF}
        \to -iK_+\left(F_{\rm reg} + iF_{\rm irr}\right). 
  \f]
 
  For a complex source, phi_P is two real standing-wave solves
  (@ref solveContinuumMixedState).

  \par Iteration
  The TDHF map is linear in the corrections, and damped iteration diverges
  near resonances (autoionising resonances, omega ~ en_b - en_a). The
  self-consistency is therefore driven by Anderson mixing
  (@ref anderson_coefficients) of the whole state, all corrections together
  with the channel amplitudes K_+ (linear in the state), which converges
  wherever the linear problem is non-singular.

  \par Usage
  As TDHF: @ref solve_core (omega), then
  - @ref A_phys: the ionisation amplitude of each open channel,
    \f$ A = K_+^* / \pi \f$, in @ref channel_list order. |A| replaces
    \f$ |\redmatel{\en\kappa}{t}{a}| \f$ in cross-sections (energy-normalised
    continuum states); the conjugate is the conventional normalisation of
    the final state to incoming waves. The phase is relative to the regular
    solution of the local channel potential: the relative phase of two
    operators in the same channel is physical, absolute phases and phases
    between channels are not.
  - @ref dV_complex (Fa, Fb): the complex \f$ \redmatel{a}{\delta V}{b} \f$;
    equals TDHF::dV below every threshold.
*/
class TDHFcntm : public TDHF {

public:
  //! Constructs for operator h_plus (with optional h_minus); see
  //! @ref TDHF::TDHF.
  TDHFcntm(const DiracOperator::TensorOperator *const h_plus,
           const HF::HartreeFock *const hf,
           const DiracOperator::TensorOperator *const h_minus = nullptr);

  //! Open channel: core-orbital index, ejected-electron kappa, and its
  //! energy en = en_a + omega
  struct Channel {
    std::size_t i_core;
    int kappa;
    double en;
  };

  /*!
    @brief Solves the (complex) TDHF equations self-consistently at omega.
    @param omega    Frequency (atomic units); its magnitude is used.
    @param max_its  Maximum number of iterations; 1 gives the first-order
                    (outgoing-wave) correction.
    @param print    If true, write convergence progress to screen.
    @details Re-solving at the same omega warm-starts from the previous
    solution; the continuum data of the open channels is rebuilt only when
    omega changes.
  */
  void solve_core(double omega, int max_its = 100, bool print = true) override;

  //! Clears the corrections and the continuum channel data
  void clear() override;

  //! Open channels at the omega of the last solve_core(), in @ref A_phys
  //! order; empty before solve_core(), or if no channel is open
  std::vector<Channel> channel_list() const;

  //! Ionisation amplitude A = K_+^*/pi of every open channel (see class
  //! description), in @ref channel_list order. Requires solve_core().
  std::vector<std::complex<double>> A_phys() const;

  //! A (see @ref A_phys) of the channel of core orbital Fa with
  //! ejected-electron kappa; zero if that channel is closed
  std::complex<double> A_phys(const DiracSpinor &Fa, int kappa) const;

  /*!
    @brief Reduced matrix element of the (complex) induced potential,
    \f$ \redmatel{a}{\delta V}{b} \f$, or the conjugate
    \f$ \redmatel{a}{\delta V^\dagger}{b} \f$ if en_b > en_a; as
    @ref TDHF::dV, complex.
    @details Real, and equal to TDHF, below every threshold. For a continuum
    bra (en_a > 0), Fa must be the energy-normalised continuum state of hole
    Fb in the V^{N-1} potential (direct y^0_bb and one-electron exchange of
    Fb removed); the source then includes the hole term V^b_0 phi of the
    rearranged channel equation, so that <Fa||t||b> + dV_complex(Fa, Fb) is
    the outgoing-wave amplitude, |A| of @ref A_phys in the phase reference
    of Fa.
  */
  std::complex<double> dV_complex(const DiracSpinor &Fa,
                                  const DiracSpinor &Fb) const;

  //! Real part of @ref dV_complex, for bound Fa and Fb (exact below every
  //! threshold, where the class acts as TDHF). For a continuum state use
  //! dV_complex.
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb) const override;
  using TDHF::dV;

private:
  // Imaginary parts of the corrections, indexed as m_X / m_Y (the real
  // parts). The Y sets store conj(Y) (see class description).
  std::vector<std::vector<DiracSpinor>> m_Xi{};
  std::vector<std::vector<DiracSpinor>> m_Yi{};

  // Continuum data of one (core orbital x channel) at fixed omega: whether
  // the channel is open (en_+ > 0 and resolvable on the grid); the
  // homogeneous pair Freg/Firr at en_+ in the local channel potential
  // (vlocal - y^0_aa + U_KS); the exchange-dressed regular solution
  // F^HF ~ Freg + K_ex Firr; and the outgoing amplitude K_+ of the latest
  // solve. Indexed as m_X.
  struct ContinuumChannel {
    bool open{false};
    DiracSpinor Freg;
    DiracSpinor Firr;
    DiracSpinor Fhf;
    double K_ex{0.0};
    std::complex<double> K{0.0};
  };
  std::vector<std::vector<ContinuumChannel>> m_channels{};
  // The omega m_channels was built at (rebuilt when it changes)
  double m_omega{-1.0};

  // (Re)shapes the imaginary sets to m_X / m_Y, zeroed
  void zero_imaginary_sets();

  // Builds m_channels for omega: for each open channel the pair Freg/Firr,
  // then F^HF and K_ex
  void prepare_channels(double omega);

  // (i_core, i_channel) of the channel of core orbital Fb with kappa
  std::pair<std::size_t, std::size_t> channel_index(const DiracSpinor &Fb,
                                                    int kappa) const;

  // One (undamped) application of the complex TDHF map to the stored
  // corrections, all (core orbital x channel x X/Y) solves task-flattened;
  // returns the convergence measure {eps, worst channel}
  std::pair<double, std::string> tdhf_core_it_complex(double omega);

  // Outgoing-wave solve of one open X channel: real and imaginary parts of
  // the previous iterate in, new iterate out; K_+ written to the channel
  void solve_channel_outgoing(DiracSpinor *X_re, DiracSpinor *X_im,
                              ContinuumChannel *channel, const DiracSpinor &Fb,
                              const DiracSpinor &hFb, const DiracSpinor &dV_re,
                              const DiracSpinor &dV_im, double omega,
                              double eps_ms) const;

  // Bound channel (closed orbital, or the Y partner of an ionised orbital),
  // as TDHF. hFb may be nullptr (no external-field source: the imaginary
  // part).
  void solve_channel_bound(DiracSpinor *dF_beta, const DiracSpinor &Fb,
                           const DiracSpinor *hFb, const DiracSpinor &dV_src,
                           double omega, dPsiType XorY, double eps_ms) const;

  // The hole term V^a_0 chi = (y^0_aa + X_a) chi of ionised orbital Fa: the
  // one-electron self-interaction moved from the Hamiltonian of the open
  // channels into their source
  DiracSpinor hole_compensation(const DiracSpinor &Fa,
                                const DiracSpinor &chi) const;

  // Convergence measure of one iteration: bound channels, the relative L2
  // change of (X_re, X_im) over the norm of ALL X channels (a negligible
  // bound remnant must not gate the iteration); open channels,
  // |dK_+|^2 / max(|K_+|^2, floor). K_old indexed as m_channels.
  std::pair<double, std::string> eps_complex(
    const std::vector<std::vector<DiracSpinor>> &Xs_re,
    const std::vector<std::vector<DiracSpinor>> &Xs_im,
    const std::vector<std::vector<std::complex<double>>> &K_old) const;

  // The full state (all corrections, real and imaginary, and the channel
  // amplitudes) as one real vector, and back: the Anderson iterate
  std::vector<double> state_vector() const;
  void set_state(const std::vector<double> &state);

public:
  TDHFcntm &operator=(const TDHFcntm &) = delete;
  TDHFcntm(const TDHFcntm &) = default;
  ~TDHFcntm() = default;
};

} // namespace ExternalField
