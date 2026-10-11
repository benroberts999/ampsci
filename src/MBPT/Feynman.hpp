#pragma once
#include "HF/HartreeFock.hpp"
#include "MBPT/RadialMatrix.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <vector>

namespace MBPT {

// FeynmanOptions class

/*!
  @brief Options for including the hole-particle interaction.

  @details
  - exclude    : not included
  - include    : included, all k
  - include_k0 : included, k=0 term only
*/
enum class HoleParticle { exclude, include, include_k0 };

//! Options for including Screening
enum class Screening { exclude, include };

//! Options for the Feynman method (see Feynman)
struct FeynmanOptions {
  //! Include all-order screening of the Coulomb interaction
  Screening screening{Screening::exclude};
  //! Include hole-particle interaction in the polarisation operator
  HoleParticle hole_particle{HoleParticle::exclude};
  //! Maximum l for internal lines (at least the maximum l in the core)
  int max_l_internal{6};
  //! Real part of the frequency, for integration along Im(w) (see best_omre)
  double omre{-0.3};
  /*!
    @brief First non-zero point of the frequency grid along Im(w).

    @details
    u=0 is always included. The [0, w0] panel is a trapezoid: w0 must be
    small compared to the pole distances (~|omre|) so the integrand is flat
    across the panel.
  */
  double w0{0.01};
  //! Ratio for the logarithmic frequency grid along Im(w)
  double w_ratio{1.5};
  /*!
    @brief Method for the Green's function at complex energy.

    @details
    - false : solve at Re(en), extend to complex energy via Dyson (resolvent)
    - true  : solve the Dirac equation directly at complex energy (Dyson still
              used at large Im(en); see Feynman::green)
  */
  bool complex_green{true};
};

/*!
  @brief Which states to include in the Green's function.

  @details
  - both    : core and excited (full Hartree-Fock Green's function)
  - core    : core states only
  - excited : excited states only (orthogonalised to the core)
*/
enum class GreenStates { both, core, excited };

//------------------------------------------------------------------------------
/*!
  @brief Direct correlation potential diagrams via the Feynman technique.

  @details
  Constructs Hartree-Fock Green's functions, the polarisation operator, and
  the direct Sigma diagrams, as matrices on a radial sub-grid. The direct
  Sigma is formed by integrating, along Im(w) at fixed Re(w) = omre, the
  product of the Green's function \f$ G(\varepsilon_v + \omega) \f$ and the
  Coulomb-polarisation-Coulomb operator, \f$ Q\Pi Q(\omega) \f$.

  Allows all-order inclusion of:
  - screening of the Coulomb interaction:
    \f$ Q\Pi Q \to Q\Pi X Q \f$, with \f$ X = [1 + iQ\Pi]^{-1} \f$
  - the hole-particle interaction, in the polarisation operator

  Exchange diagrams are not formed here (see Goldstone).
*/
class Feynman {

  // Pointer to HF potential, and HF core
  const HF::HartreeFock *m_HF;
  // Pointer to shared radial grid (full grid)
  std::shared_ptr<const Grid> m_grid;
  // Parameters of the sub-grid: initial/final points, stride
  std::size_t m_i0, m_stride, m_subgrid_points;

  // maximum kappa_index appearing in the core
  std::size_t m_max_ki_core;
  // maximum kappa index to include for internal lines
  std::size_t m_max_ki;
  // maximum multipolarity, k
  int m_max_k;
  // Lowest n to polarise in polarisation operator
  int m_min_core_n;
  // Include the small (g) spinor components: in the internal lines (Green's
  // functions, exchange and projection matrices, polarisation loop) and in
  // Sigma itself. Changes Sigma by ~0.05%, but costs ~8x in the linear
  // algebra. Always on with Breit (Vx then mixes f and g)
  bool m_include_G;
  // real part of frequency for integration
  double m_omre;
  // Parameters of the log grid along Im(w): first point, ratio
  double m_w0;
  double m_wratio;
  // Frequency quadrature along Im(w): points u_i (u_0 = 0 explicitly) and
  // weights W_i in linear measure: int_0^inf F(u) du =~ sum_i W_i F(u_i).
  // Weights include the [0, w0] end panel and the u^-3 large-u tail.
  std::vector<double> m_wgrid_points{};
  std::vector<double> m_wgrid_weights{};

  // Option: include hole-particle interaction into polarisation operator
  bool m_hole_particle;
  // Option to include higher-k in hp interaction
  bool m_include_higher_order_hp;
  // Include screening correction to Coulomb line in Q*Pi*Q => Q*Pi*X*Q
  bool m_screen_Coulomb;

  // Coulomb operator; include right integration measure
  // q_ij := (r>^k/r<^{k+1}) * dr_j
  std::vector<ComplexRMatrix> m_qk{};

  // Core projection operators (one for each core state)
  std::vector<ComplexGMatrix> m_pa{};
  // Summed core projection operator, P = sum_a |a><a| (one for each kappa)
  std::vector<ComplexGMatrix> m_Pcore{};
  // Hartree-Fock Exchange matrix (one for each kappa)
  std::vector<GMatrix> m_Vx_kappa{};

  // Effective spinless Q*Pi*Q operator: for each imaginary omega, and each k
  LinAlg::Matrix<ComplexRMatrix> m_qpiq_wk{};

  // Method for complex-energy Green's function (see FeynmanOptions)
  // false: Dyson method (solve at real energy, correct to complex);
  // true: solve Dirac equation directly at complex energy
  bool m_Complex_green_method;

public:
  /*!
    @brief Sets up the Feynman method on given sub-grid, and forms Q*Pi*Q.

    @details
    Forms the radial Coulomb (q^k), core projection, and Hartree-Fock
    exchange matrices, and the frequency quadrature along Im(w). Q*Pi*Q is
    the expensive step: only needed for Sigma_direct (the Green's functions
    etc. do not need it), and can be skipped (@p form_qpq = false) and formed
    later with form_qpiq().

    @param vHF         Hartree-Fock (core, potential, grid). Pointer is
                       stored: must outlive this object
    @param i0          Sub-grid: index of first point on the full grid
    @param stride      Sub-grid stride: every stride-th point of full grid
    @param size        Number of sub-grid points
    @param options     Screening, hole-particle, internal lines, frequency
                       grid, Green's function method (see FeynmanOptions)
    @param n_min_core  Minimum n of core states in the polarisation loop
    @param include_G   Include small (g) components throughout. Always on if
                       @p vHF includes Breit
    @param verbose     Print summary of options
    @param ident       File prefix for the Q*Pi*Q disk cache ("" or "false":
                       no cache)
    @param form_qpq    Form Q*Pi*Q on construction
  */
  Feynman(const HF::HartreeFock *vHF, std::size_t i0, std::size_t stride,
          std::size_t size, const FeynmanOptions &options, int n_min_core,
          bool include_G, bool verbose = true, const std::string &ident = "",
          bool form_qpq = true);

  //! True if Coulomb screening is included
  bool screening() const { return m_screen_Coulomb; }
  //! True if hole-particle interaction is included
  bool hole_particle() const { return m_hole_particle; }

  //! Returns stride used for sub-grid
  std::size_t stride() const { return m_stride; }

  //! Minimum n of core states in the polarisation loop
  int n_min() const { return m_min_core_n; }
  //! Maximum multipolarity, k
  int max_k() const { return m_max_k; }
  //! Real part of frequency used for integration
  double omre() const { return m_omre; }
  //! First non-zero point of the frequency grid along Im(w)
  double w0() const { return m_w0; }
  //! Ratio of the logarithmic frequency grid along Im(w)
  double wratio() const { return m_wratio; }
  //! Maximum l for internal lines
  int lmax() const { return Angular::kindex_to_l(m_max_ki); }

  /*!
    @brief Green's function for given kappa, at complex energy @p en.

    @details
    For @p states = both, the full Hartree-Fock Green's function. With
    complex_green set (see FeynmanOptions): solved directly at complex energy
    for small Im(en), and by the Dyson method at large Im(en) (where it is
    consistent with the discrete radial sums). Otherwise, Dyson method
    throughout.

    @param kappa   Dirac quantum number
    @param en      Complex energy
    @param states  Which states to include:
                   - both    : core and excited
                   - core    : core only, sum over core states a of
                               |a><a|/(en - e_a)
                   - excited : excited only (orthogonalised to the core)
  */
  ComplexGMatrix green(int kappa, std::complex<double> en,
                       GreenStates states = GreenStates::both) const;

  /*!
    @brief Green's function, solving the Dirac equation directly at complex en.

    @details
    cf. green(), which mixes this with the Dyson method.
  */
  ComplexGMatrix green_complex_dirac(int kappa, std::complex<double> en) const;

  /*!
    @brief Green's function by the Dyson (resolvent) method.

    @details
    Solves at Re(en), then shifts to complex energy. cf. green(), which mixes
    the two methods.
  */
  ComplexGMatrix green_dyson(int kappa, std::complex<double> en) const;

  /*!
    @brief Green's function for given kappa, as a sum over an explicit basis.

    @details
    \f$ G = \sum_n \ket{n}\bra{n}/(\varepsilon - \varepsilon_n) \f$, over
    states in @p basis with given kappa.
  */
  ComplexGMatrix green_basis(int kappa, std::complex<double> en,
                             const std::vector<DiracSpinor> &basis) const;

  /*!
    @brief Polarisation operator pi^k(w), for each k = 0..max_k separately.

    @details
    Not summed over k:
    \f[
      \Pi^k(\omega) = i \sum_{a,n} \frac{(C^k_{an})^2}{2k+1} \,
      F_a^\dagger \left[ G^{\rm ex}_n(\varepsilon_a - \omega)
      + G^{\rm ex}_n(\varepsilon_a + \omega) \right] F_a,
    \f]
    with core states a (n >= n_min) and excited-state Green's functions for
    each internal kappa_n. Green's functions are computed once and re-used
    across all k.

    @param omega          Complex frequency
    @param hole_particle  Include the hole-particle interaction (independent
                          of the option set on construction)
  */
  std::vector<ComplexRMatrix> polarisation_each_k(std::complex<double> omega,
                                                  bool hole_particle) const;

  /*!
    @brief Polarisation operator pi^k(w) on the frequency grid, for each k.

    @details
    This is the expensive stage of forming Q*Pi*Q. Nothing is stored: pass
    the result to form_qpiq(pi_wk) to form Q*Pi*Q (for this object, or for
    another that differs only in screening). Depends on the hole-particle
    option, not on screening.

    @return pi^k(w) at every point of the frequency grid, indexed [iw][k].
  */
  std::vector<std::vector<ComplexRMatrix>> polarisation_wk() const;

  /*!
    @brief Forms Q*Pi*Q along the frequency grid (required by Sigma_direct).

    @details
    With @p ident, reads from the disk cache if a matching file exists, else
    forms and writes it (see read_qpiq/write_qpiq).
  */
  void form_qpiq(const std::string &ident = "");

  /*!
    @brief Forms Q*Pi*Q (with screening, if set) from given polarisation op.

    @details
    @p pi_wk is from polarisation_wk(), and may be shared between Feynman
    objects that differ only in screening. With screening,
    \f$ Q\Pi Q \to Q\Pi X Q \f$, with \f$ X = [1 + iQ\Pi]^{-1} \f$.

    @warning @p pi_wk must be on the same sub-grid and frequency grid
             (checked only by assert).
  */
  void form_qpiq(const std::vector<std::vector<ComplexRMatrix>> &pi_wk);

  //! True once Q*Pi*Q has been formed
  bool has_qpiq() const { return m_qpiq_wk.size() != 0; }

  /*!
    @brief Reads Q*Pi*Q from the disk cache for @p ident (file prefix).

    @return False if no matching file (or @p ident is "" or "false").
  */
  bool read_qpiq(const std::string &ident);

  //! Writes Q*Pi*Q to the disk cache for ident; false if not written
  bool write_qpiq(const std::string &ident);

  /*!
    @brief Direct part of the correlation potential, for kappa_v at en_v.

    @details
    Summed over multipoles k (unless @p k given) and internal kappa. See
    class description.

    @param kappa_v  Valence Dirac quantum number
    @param en_v     Energy at which Sigma is formed
    @param k        If given, only this multipole (zero matrix if out of
                    range)

    @warning Q*Pi*Q must be formed first (see form_qpiq); checked only by
             assert.
  */
  GMatrix Sigma_direct(int kappa_v, double en_v,
                       std::optional<int> k = {}) const;

  /*!
    @brief Direct part of correlation potential, for each multipole k.

    @details
    Not summed: Sigma_d = sum of these. Same cost as a single Sigma_direct
    call, since the Green's functions are shared by all k.

    @warning Q*Pi*Q must be formed first (see form_qpiq).
  */
  std::vector<GMatrix> Sigma_direct_each_k(int kappa_v, double en_v) const;

  /*!
    @brief Returns (reference to) radial Coulomb matrix, q^k.

    @details
    \f$ q^k_{ij} = (r_<^k / r_>^{k+1}) \, dr_j \f$: includes the right-hand
    integration measure (drj) only.
  */
  const ComplexRMatrix &get_qk(int k) const { return m_qk.at(std::size_t(k)); }

  //! Returns (ref to) radial exchange matrix Vx_kappa. Nb: includes dri*drj
  const GMatrix &get_Vx_kappa(int kappa) const {
    return m_Vx_kappa.at(Angular::kappa_to_kindex(kappa));
  }

  /*!
    @brief Returns (ref to) summed core projection operator, for given kappa.

    @details
    \f$ P = \sum_a \ket{a}\bra{a} \f$. No integration measure.
  */
  const ComplexGMatrix &get_Pcore(int kappa) const {
    return m_Pcore.at(Angular::kappa_to_kindex(kappa));
  }

private:
  // forms Qk matrices, as well as dri, drj
  void form_qk();
  // Forms core projection operators
  void form_pa();
  // Forms HF exchange potential matrix
  void form_vx();
  // Sets up the frequency quadrature (points + weights) along Im(w):
  // u = 0 explicitly (the integrand is finite, generically maximal, there),
  // composite Simpson's rule in ln(u) on the log grid [w0, wmax], trapezoid
  // panel for [0, w0], and u^-3 tail correction beyond wmax
  void form_w_quadrature(double w0, double wratio);
  // Disk-cache filename for Q*Pi*Q (encodes the options it depends on);
  // "" if caching is disabled for this ident
  std::string qpiq_filename(const std::string &ident) const;

  bool readwrite_qpiq(IO::FRW::RoW rw, const std::string &fname);

  // Screening factor X = [1 + i qk*pik]^-1
  ComplexRMatrix X_screen(const ComplexRMatrix &pik,
                          const ComplexRMatrix &qdri) const;

  // Forms single "Green's function" contribution f*|ket><bra|
  // (f is usually 1/(e-en))
  ComplexGMatrix green_single(const DiracSpinor &ket, const DiracSpinor &bra,
                              const std::complex<double> f) const;

  // Forms full Hartree-Fock green's function. If Fc_hp is given,
  // includes hole-particle contribution (Fc is polarised core state).
  // Uses Dyson method to extend to complex energies.
  ComplexGMatrix green_hf(int kappa, std::complex<double> en,
                          const DiracSpinor *Fc_hp = nullptr) const;

  // Forms full Hartree-Fock green's function. If Fc_hp is given,
  // includes hole-particle contribution (Fc is polarised core state).
  // Uses "Complex Dirac" method for complex energies.
  ComplexGMatrix
  green_hf_complex_dirac(int kappa, std::complex<double> en,
                         const DiracSpinor *Fc_hp = nullptr) const;

  // Given Gr(wr), and wi (wr is real), returns Gr(wr + i*wi)
  ComplexGMatrix green_to_complex(const ComplexGMatrix &Gr,
                                  double om_imag) const;

  // Green's function, only due to core states
  ComplexGMatrix green_core(int kappa, std::complex<double> en) const;

  // Greens function, only due to excited states (by orthog to core)
  ComplexGMatrix green_excited(int kappa, std::complex<double> en,
                               const DiracSpinor *Fc_hp = nullptr) const;

  // Projects core states out of Green's fn from both sides:
  // G -> (1-P) G (1-P), with P = sum_a |a><a| dr (discrete brakets)
  ComplexGMatrix orthogonalise_wrt_core(const ComplexGMatrix &g_in,
                                        int kappa) const;

  // Replaces each sampled cell of a complex-energy Green matrix by its
  // average over the dr*dr cell (near-diagonal cells are otherwise
  // unresolved on the sub-grid); factors from the solutions' local decay
  void cell_average(ComplexGMatrix *G, const DiracSpinor &x0,
                    const DiracSpinor &Ix0, const DiracSpinor &xI,
                    const DiracSpinor &IxI) const;

  // Given Dirac solutions regular at 0 (x0) and infinity (xI), forms "local"
  // Green's function, with all four spinor components
  GMatrix construct_green_g0(const DiracSpinor &x0, const DiracSpinor &xI,
                             const double w) const;

  // Given Dirac solutions regular at 0 (x0 + i*Ix0) and infinity (xI + i*IxI),
  // forms "local" Green's function (complex).
  ComplexGMatrix construct_green_g0(const DiracSpinor &x0,
                                    const DiracSpinor &Ix0,
                                    const DiracSpinor &xI,
                                    const DiracSpinor &IxI,
                                    const std::complex<double> w) const;

  // Calculates hole-particle potential (1-P)Vhp(1-P)
  GMatrix calculate_Vhp(int kappa, const DiracSpinor &Fc) const;

  // [[nodiscard]] const ComplexGMatrix &get_dri() const { return m_dri; }
  // [[nodiscard]] const ComplexGMatrix &get_drj() const { return m_drj; }

public:
  Feynman &operator=(const Feynman &) = default;
  Feynman(const Feynman &) = default;
  ~Feynman() = default;
};

/*!
  @brief Best real part of the frequency, omre, for the direct diagram.

  @details
  The value furthest from any pole of the u=0 integrand. With
  \f$ \Delta = \varepsilon_{\rm exc} - \varepsilon_{\rm core} \f$
  (lowest excited minus highest core energy), the window
  \f$ (-\Delta, 0) \f$ is free of true poles; inside it lie only the weak
  fictitious poles at core-core energy differences (imperfect core
  subtraction in the polarisation-loop Gex) and, for non-lowest valence
  states, the small valence-valence differences (valence-line G). Only the
  lowest few physical excited states enter, so the valence list suffices: no
  basis/spectrum needed (a typical valence energy is assumed if the list is
  empty). One omre serves all valence states.

  @param core     Core states
  @param valence  Valence states (lowest sets the window)
  @param print    Print Delta, the in-window poles, and the chosen omre

  @return Midpoint of the largest pole-free gap.
*/
double best_omre(const std::vector<DiracSpinor> &core,
                 const std::vector<DiracSpinor> &valence, bool print = false);

} // namespace MBPT