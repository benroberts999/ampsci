#pragma once
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/Goldstone.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <optional>
#include <vector>

namespace MBPT {

//! Stored Sigma matrix, with the parameters it was formed with
struct SigmaData {
  //! Dirac quantum number of the state Sigma was formed for
  int kappa;
  //! Energy Sigma was formed at
  double en;
  //! Sigma matrix, on the sub-grid
  SpinorMatrix<double> Sigma;
  //! Principal quantum number of the state Sigma was formed for (0: none given)
  int n{0};
  //! Scaling factor, Sigma -> lambda*Sigma; not applied to the stored matrix
  double lambda{1.0};
  //! fk screening factors used for this Sigma (empty if none)
  std::vector<double> fk{};
  //! Coefficient of large-r tail, Sigma = a_eff/r^4 beyond sub-grid (0: none)
  double a_eff{0.0};
};

/*!
  @brief Fits the coefficient a_eff of the large-r polarisation tail of Sigma.

  @details
  Large-r tail: \f$ \Sigma(r) \to a_{\rm eff}/r^4 \f$.

  Averages \f$ r^4 (\Sigma F_v)(r)/F_v(r) \f$ (large component) over @p n_fit
  sub-grid points: the outermost ones where
  \f$ |F_v| \geq f_{\rm cut} \max|F_v| \f$ (stopping at a node). Further out,
  Fv is exponentially small and the ratio is dominated by numerical noise in
  Sigma. In the fit region (outside the core, typically a few to ~15 aB)
  \f$ r^4 (\Sigma F_v)/F_v \f$ is already close to constant.

  @return a_eff; zero if there are no points in the fit region.

  @note Can be numerically unstable, should be treated as a test
*/
double fit_tail_coefficient(const GMatrix &Sigma, const DiracSpinor &Fv,
                            std::size_t n_fit = 10, double f_cut = 1.0e-2);

/*!
  @brief Adds the large-r polarisation tail of Sigma to @p SigmaFn.

  @details
  Adds \f$ (a_{\rm eff}/r^4) F_n \f$ for r beyond the Sigma sub-grid (where
  \f$ \Sigma F_n \f$ is zero). Treated as a local potential. Acts on the small
  component only if @p Sigma includes g parts. Does nothing if @p a_eff is
  zero.
*/
void add_tail(DiracSpinor *SigmaFn, const GMatrix &Sigma, const DiracSpinor &Fn,
              double a_eff);

/*!
  @brief Method used to form the direct Sigma diagrams.

  @details
  Goldstone (basis sum) or Feynman (Green's function). Exchange diagrams
  always use Goldstone.
*/
enum class SigmaMethod { Goldstone, Feynman };

//! Sigma sub-grid parameters: range (aB) and stride (of full grid)
struct rgrid_params {
  double r0{1.0e-4};
  double rmax{30.0};
  std::size_t stride{4};
};

//==============================================================================

/*!
  @brief Correlation potential, Sigma, for valence states.

  @details
  Core-valence correlations, stored as a non-local matrix for each kappa, on
  a radial sub-grid (every stride-th point of the full grid, between r0 and
  rmax). Formed at a given energy (see formSigma). Read from file if it
  exists, so it need only be formed once (see write).

  Direct diagrams are formed using either the Goldstone (basis sum) or Feynman
  (Green's function) technique; the Feynman technique allows all-order
  inclusion of screening and the hole-particle interaction. Exchange diagrams
  always use Goldstone, with optional effective screening factors fk.

  The action of Sigma on a state (SigmaFv) includes:
  - lambda scaling, \f$ \Sigma \to \lambda\Sigma \f$ (see scale_Sigma)
  - ladder correction, \f$ \Sigma_L \f$, if a ladder file was given
  - large-r tail, \f$ a_{\rm eff}/r^4 \f$, beyond the sub-grid, if
    tail_extrapolation is set (experimental; see fit_tail_coefficient)

  Optionally also forms the energy derivative,
  \f$ \partial\Sigma/\partial\varepsilon \f$ (see dSigmaFv).

  @note Not thread safe: formSigma modifies internal state.
*/
class CorrelationPotential {
  const HF::HartreeFock *m_HF;
  std::vector<DiracSpinor> m_basis; // so we can delay goldstone construction
  std::vector<SigmaData> m_Sigmas{};
  double m_r0, m_rmax;
  std::size_t m_stride;
  std::size_t m_i0, m_size; // need?

  SigmaMethod m_method;
  int m_n_min_core;
  // int m_n_min_core_F;
  bool m_includeG;
  bool m_includeBreit_b2;
  int m_n_max_breit;

  std::optional<Goldstone> m_Gold{};

  FeynmanOptions m_Foptions;
  bool m_calculate_fk; // if not, need fk and etak
  std::vector<double> m_fk;
  std::vector<double> m_etak;
  // Apply fk to both Coulomb lines of the exchange diagrams (else outer only)
  bool m_fk_both_lines;

  std::optional<Feynman> m_Fy{};

  // These are only for calculating fk and eta
  std::optional<Feynman> m_Fy0{};
  std::optional<Feynman> m_FyX{};
  std::optional<Feynman> m_FyH{};

  std::string m_fname{};

  // Ladder correction, Sigma_L: read from file (produced by the Ladder{}
  // block/driver - see MBPT::ladder). Stored separately from the base Sigma
  // (may have its own sub-grid and include_G); added in SigmaFv, and scaled
  // by the same lambda as the base Sigma.
  std::vector<SigmaData> m_Sigma_L{};
  std::string m_ladder_file{};

  // Energy derivative, dSigma/dE (forward difference; see m_delta_en).
  // Formed alongside Sigma when m_form_derivative is set; appended to the
  // sigma file (older files simply have none)
  bool m_form_derivative{false};
  std::vector<SigmaData> m_dSigma{};
  static constexpr double m_delta_en = 0.01;

  // Experimental: extrapolate Sigma beyond the sub-grid as a_eff/r^4, with
  // a_eff fitted once for each Sigma from the valence state it was formed for
  // (see fit_tail_coefficient). Sigmas formed without a valence state get no
  // tail.
  bool m_tail_extrapolation{false};

public:
  /*!
    @brief Reads Sigma from file if it exists; otherwise sets up to form it.

    @details
    If @p fname exists, Sigma is read from it. Otherwise, sets up the
    Goldstone/Feynman method; Sigma itself is only formed by formSigma. The
    ladder file, if given, is read in either case.

    @param fname               Sigma file (read if it exists; see write)
    @param vHF                 Hartree-Fock (core, grid, Breit). Pointer is
                               stored: must outlive this object
    @param basis               Basis for the Goldstone sums (and print_de)
    @param r0                  Sub-grid first point (aB)
    @param rmax                Sub-grid last point (aB)
    @param stride              Sub-grid stride: every stride-th point of the
                               full grid
    @param n_min_core          Minimum n of core states included in
                               polarisation (lower core states are frozen)
    @param method              Goldstone or Feynman, for direct diagrams
    @param include_g           Include small-component (fg, gf, gg) parts of
                               Sigma
    @param include_Breit_b2    Include two-body Breit (B2) correction (only if
                               @p vHF includes Breit)
    @param n_max_breit         Maximum n of excited states in the B2
                               correction (0: all)
    @param Foptions            Options for the Feynman method (screening,
                               hole-particle, etc.)
    @param calculate_fk        Calculate fk screening factors for exchange,
                               for each state (Feynman, with screening); else
                               @p fk used
    @param fk                  Effective screening factors, for each
                               multipole k
    @param etak                Effective hole-particle factors, for each k
                               (Goldstone only)
    @param ladder_file         File with ladder correction, Sigma_L (empty:
                               none)
    @param form_derivative     Also form dSigma/dE (see dSigmaFv)
    @param fk_both_lines       Apply fk to both Coulomb lines of the exchange
                               diagrams (else outer only)
    @param tail_extrapolation  Extrapolate Sigma beyond the sub-grid as
                               a_eff/r^4 (experimental; see
                               fit_tail_coefficient)

    @note If read from file, the sub-grid and @p include_g are taken from the
          file, not the inputs.
  */
  CorrelationPotential(
    const std::string &fname, const HF::HartreeFock *vHF,
    const std::vector<DiracSpinor> &basis, double r0, double rmax,
    std::size_t stride, int n_min_core, SigmaMethod method,
    bool include_g = false, bool include_Breit_b2 = false, int n_max_breit = 0,
    const FeynmanOptions &Foptions = {}, bool calculate_fk = true,
    const std::vector<double> &fk = {}, const std::vector<double> &etak = {},
    const std::string &ladder_file = "", bool form_derivative = false,
    bool fk_both_lines = false, bool tail_extrapolation = false);

  /*!
    @brief Forms and stores Sigma for given kappa, at energy @p ev.

    @details
    If Sigma already exists for @p kappa and @p n (any n if @p n <= 0), it is
    not re-formed; only the tail fit and dSigma/dE are added if required and
    missing (e.g., Sigma read from file).

    @param kappa  Dirac quantum number
    @param ev     Energy at which Sigma is formed
    @param n      Principal quantum number: allows distinct Sigma for each n
                  of the same kappa
    @param Fv     State Sigma is formed for (optional): used to calculate fk
                  screening factors, fit the large-r tail, and print the
                  energy shift. Must have kappa = @p kappa.

    @warning Not thread safe.
  */
  void formSigma(int kappa, double ev, int n, const DiracSpinor *Fv = nullptr);

  //! True if no Sigma matrices are stored
  bool empty() const { return m_Sigmas.empty(); }

  /*!
    @brief Pointer to the stored Sigma matrix for given kappa (and n).

    @details
    See get() for the choice of n.

    @return Pointer to the Sigma matrix; nullptr if not present.

    @note The raw matrix: does not include lambda, the ladder part, or the
          large-r tail (use SigmaFv for those).
  */
  const GMatrix *getSigma(int kappa, int n = 0) const;

  //! Scaling factor, lambda, for given kappa (and n); 1 if not present
  double getLambda(int kappa, int n = 0) const;

  //! Removes all stored Sigma matrices (not ladder or dSigma/dE)
  void clear() { m_Sigmas.clear(); }

  /*!
    @brief Returns Sigma|Fv>, including lambda scaling, ladder, and tail.

    @details
    \f$ \lambda (\Sigma + \Sigma_L) \ket{F_v} \f$: scaled by lambda, and
    includes the ladder correction \f$ \Sigma_L \f$ if present. Includes the
    large-r tail if it was fitted (see tail_extrapolation option).

    @return Sigma|Fv>; |0> if neither Sigma nor Sigma_L exists for kappa_v.
  */
  DiracSpinor SigmaFv(const DiracSpinor &Fv) const;
  //! Same as SigmaFv
  DiracSpinor operator()(const DiracSpinor &Fv) const { return SigmaFv(Fv); }

  /*!
    @brief Prints MBPT(2) energy correction for each valence state.

    @details
    Direct, exchange, and total, and \f$ \langle v|\Sigma|v \rangle \f$ using
    the stored Sigma matrix.
  */
  void print_de(const std::vector<DiracSpinor> &valence);

  /*!
    @brief Prints a_eff for each valence state with a fitted large-r tail.

    @details
    Also prints \f$ \langle v|\Sigma|v \rangle \f$ without (S0) and with (S)
    the tail. Prints nothing if tail_extrapolation is off.
  */
  void print_tail(const std::vector<DiracSpinor> &valence) const;

  //! True if any dSigma/dE matrices are present (see form_derivative)
  bool has_derivative() const { return !m_dSigma.empty(); }

  //! Pointer to dSigma/dE data for given kappa (and n); nullptr if not present
  const SigmaData *get_derivative(int kappa, int n = 0) const;

  /*!
    @brief Returns lambda * dSigma/dE |Fv>.

    @details
    \f$ \lambda \, (\partial\Sigma/\partial\varepsilon) \ket{F_v} \f$.
    The energy derivative of Sigma, formed by finite difference alongside
    Sigma itself (option form_derivative), so it corresponds to the actual
    method used (Goldstone/Feynman, screening, etc.). Scaled by the same
    lambda as the base Sigma.

    @return lambda * dSigma/dE |Fv>; |0> if no derivative exists for this
            kappa.

    @note The ladder correction (Sigma_L) has no energy derivative: it is
          included in SigmaFv() but not here.
  */
  DiracSpinor dSigmaFv(const DiracSpinor &Fv) const;

  /*!
    @brief Stores scaling factors, lambda, for each stored Sigma.

    @details
    \f$ \Sigma \to \lambda\Sigma \f$. One factor for each stored Sigma, in the
    order they were formed.
  */
  void scale_Sigma(const std::vector<double> &lambdas);

  /*!
    @brief Stores scaling factor, lambda, for Sigma with given kappa and n.

    @details
    \f$ \Sigma \to \lambda\Sigma \f$. If @p n <= 0, for all Sigma with given
    kappa.
  */
  void scale_Sigma(double lambda, int kappa, int n = 0);

  //! Prints the scaling factors to screen
  void print_scaling() const;

  //! Prints kappa, energy, lambda, and if dSigma/dE exists, for each Sigma
  void print_info() const;

  /*!
    @brief Short string identifying the method (main differences only).

    @details
    E.g., "Feynman, all-order" or "Goldstone". Intended for file-cache keys.
  */
  std::string method_string() const;

  /*!
    @brief Fitting (scaling) factors as a string, "kappa=value,".

    @details
    Rounded to 4 dp; empty if none scaled. For file-cache keys (e.g., the
    ci-file hash).
  */
  std::string lambda_string() const;

  /*!
    @brief Average of the stored fk screening factors, over l up to @p l_max.

    @details
    Averaged over the lowest Sigma of each l. For re-use of the Sigma_1
    screening in Sigma_2 (CI), where there is one effective fk per Coulomb
    line but no unique state: the low-l factors are the relevant ones for the
    CI valence space (higher-l factors can differ significantly, e.g., f/g at
    k = 0, but such states rarely contribute). Weighted by l, not kappa: the
    fine-structure pair of each l is averaged first, then the mean is taken
    over the ls (so s counts the same as p, d, ...). Truncated to the shortest
    stored list.

    @param l_max  Maximum l included in the average (CI uses 2: s, p, d).

    @return Averaged fk, for each k; empty if none stored.
  */
  std::vector<double> average_fk(int l_max = 2) const;

  //! Method used for the direct diagrams: Goldstone or Feynman
  SigmaMethod method() const { return m_method; }

  /*!
    @brief Feynman object used for the direct diagrams (sub-grid, options,
           Q*Pi*Q), if formed.

    @details
    It is formed with Sigma; not when Sigma is read from file (see
    make_feynman). nullptr if not formed, or not the Feynman method.
  */
  const Feynman *feynman() const { return m_Fy ? &*m_Fy : nullptr; }

  /*!
    @brief Constructs a Feynman object with the sub-grid and parameters of
           this Sigma, with the given screening and hole-particle options.

    @details
    Inherited from this Sigma (the Correlations input): the radial sub-grid
    (rmin, rmax, stride), lmax for internal lines, n_min_core for the
    polarisation loop, include_G (small components), and the frequency
    contour (real_omega, imag_omega) and complex-energy Green's function
    method. Only screening and hole-particle are set here.

    Q*Pi*Q is read from the disk cache used when forming Sigma, if a file
    with these options is present (its name encodes them), else formed and
    written. For when feynman() is not available (Sigma read from file), or
    has other options. Empty if not the Feynman method.
  */
  std::optional<Feynman> make_feynman(Screening screening,
                                      HoleParticle hole_particle) const;

  //! Prints the sub-grid parameters to screen
  void print_subGrid() const;

  /*!
    @brief Writes Sigma matrices to file @p fname.

    @details
    Includes the sub-grid, dSigma/dE, and fk. Not written: lambda, a_eff
    (re-fitted in formSigma), or the ladder correction.
  */
  void write(const std::string &fname) { read_write(fname, IO::FRW::write); }

  /*!
    @brief Pointer to the stored Sigma data for given kappa (and n).

    @details
    If @p n <= 0, the first stored Sigma with given kappa. Otherwise, the one
    with given n; if none, the nearest lower n, else the first with given
    kappa.

    @return Pointer to Sigma data (matrix, energy formed at, lambda, etc.);
            nullptr if not present.

    @note Sigma matrix is raw (no lambda, ladder, or tail; see SigmaFv).
  */
  const SigmaData *get(int kappa, int n = 0) const;

private:
  bool read_write(const std::string &fname, IO::FRW::RoW rw);
  void setup_Feynman();
  std::vector<double> calculate_fk(double ev, const DiracSpinor &v) const;
  std::vector<double> calculate_etak(double ev, const DiracSpinor &v) const;
  const SigmaData *get_ladder(int kappa, int n = 0) const;

  // given_fk: screening factors to use (from state_fk); if nullptr,
  // calculated internally (if applicable). print = false: silent (e.g., the
  // extra evaluation for the derivative - only the used Sigma is reported)
  GMatrix formSigma_F(int kappa, double ev, const DiracSpinor *Fv = nullptr,
                      const std::vector<double> *given_fk = nullptr,
                      bool print = true);
  GMatrix formSigma_G(int kappa, double ev, const DiracSpinor *Fv = nullptr,
                      bool print = true);

  // Calculates the fk screening factors for this state, once (with print;
  // stores first set in m_fk). nullopt if not applicable (not
  // Feynman-screening, fk given manually, or no Fv)
  std::optional<std::vector<double>> state_fk(double ev, const DiracSpinor *Fv);

  // Forms dSigma/dE for given kappa by forward finite difference (step
  // m_delta_en): one extra Sigma evaluation at ev + delta, re-using the
  // base Sigma matrix and the same fk; stores in m_dSigma
  void form_derivative(int kappa, double ev, int n, const DiracSpinor *Fv,
                       const GMatrix &Sigma0,
                       const std::vector<double> *given_fk = nullptr);

  // If tail_extrapolation is set and not yet fitted, fits a_eff for this
  // Sigma from the valence state it was formed for (no tail if Fv is null)
  void fit_tail(SigmaData *Sig, const DiracSpinor *Fv) const;

public:
  CorrelationPotential &operator=(const CorrelationPotential &) = default;
  CorrelationPotential(const CorrelationPotential &) = default;
  ~CorrelationPotential() = default;

  //
};

} // namespace MBPT
