#pragma once
#include "Angular/SixJTable.hpp"
#include "HF/HartreeFock.hpp"
#include "MBPT/RadialMatrix.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <optional>
#include <vector>

namespace MBPT {

// FeynmanOptions class

//! Options for including hole-particle interaction. include mean all k;
//! include_k0 means k=0 term only
enum class HoleParticle { exclude, include, include_k0 };

//! Options for including Screening
enum class Screening { exclude, include };

struct FeynmanOptions {
  Screening screening{Screening::exclude};
  HoleParticle hole_particle{HoleParticle::exclude};
  int max_l_internal{6};
  double omre{-0.3};
  // nb: w0 = first non-zero point along Im(w) (u=0 is always included).
  // The [0, w0] panel is a trapezoid: w0 must be small compared to the
  // pole distances (~|omre|) so the integrand is flat across the panel.
  double w0{0.01};
  double w_ratio{1.5};
  // Method for Green's function at complex energy:
  // false: solve at Re(en), extend to complex energy via Dyson (resolvent);
  // true: solve the Dirac equation directly at complex energy
  bool complex_green{true};
  // int n_min_core{1};
};

//! Which states to include in Green's function:
enum class GreenStates { both, core, excited };

//------------------------------------------------------------------------------
//! Class to construct Feynman diagrams, Green's functions and polarisation op.
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
  //! Construct Feynman diagram
  //! @details rgrid_params={r0, rmax, stride}; omre is real part of frequency,
  //! w_0 is initial point along imaginary freq. axis;
  //! w_ratio is ratio used for logarithic omega grid (integration);
  //! scr_option and hp_option are screening and hole-particle interactions;
  //! max_l is maximum l to include for internal lines (Green's functions);
  //! n_min_core is minimum n to include in polarisation loop;
  //! ident is the file prefix for the Q*Pi*Q disk cache ("" or "false": no
  //! cache); form_qpq=false skips forming Q*Pi*Q (the expensive step: only
  //! needed for Sigma_direct; the Green's functions etc. do not need it),
  //! which can be formed later with form_qpiq()
  Feynman(const HF::HartreeFock *vHF, std::size_t i0, std::size_t stride,
          std::size_t size, const FeynmanOptions &options, int n_min_core,
          bool include_G, bool verbose = true, const std::string &ident = "",
          bool form_qpq = true);

  bool screening() const { return m_screen_Coulomb; }
  bool hole_particle() const { return m_hole_particle; }

  //! Returns stride used for sub-grid
  std::size_t stride() const { return m_stride; }

  int n_min() const { return m_min_core_n; }
  int max_k() const { return m_max_k; }
  double omre() const { return m_omre; }
  double w0() const { return m_w0; }
  double wratio() const { return m_wratio; }
  int lmax() const { return Angular::kindex_to_l(m_max_ki); }

  //! Calculates Green's function for kappa, and complex energy
  ComplexGMatrix green(int kappa, std::complex<double> en,
                       GreenStates states = GreenStates::both) const;

  //! Green's function by solving the Dirac equation directly at complex
  //! energy (cf. green(), which uses the Dyson method unless complex_green set)
  ComplexGMatrix green_complex_dirac(int kappa, std::complex<double> en) const;

  //! Green's function by the Dyson (resolvent) method: solve at Re(en), then
  //! shift to complex energy (cf. green(), which mixes the two methods)
  ComplexGMatrix green_dyson(int kappa, std::complex<double> en) const;

  // Uses explicit basis
  ComplexGMatrix green_basis(int kappa, std::complex<double> en,
                             const std::vector<DiracSpinor> &basis) const;

  //! Polarisation operator pi^k(w), for each k = 0..max_k separately (not
  //! summed). Green's functions are computed once and re-used across all k.
  std::vector<ComplexRMatrix> polarisation_each_k(std::complex<double> omega,
                                                  bool hole_particle) const;

  //! Calculates and returns the polarisation operator pi^k(w) at every point
  //! of the frequency grid, for each k (indexed [iw][k]). Nothing is stored:
  //! pass the result to form_qpiq(pi_wk) to form Q*Pi*Q (for this object, or
  //! for another that differs only in screening). This is the expensive stage
  //! of forming Q*Pi*Q; depends on the hole-particle option, not on screening
  std::vector<std::vector<ComplexRMatrix>> polarisation_wk() const;

  //! Forms Q*Pi*Q along the frequency grid (required by Sigma_direct).
  //! With ident, reads from the disk cache if a matching file exists, else
  //! forms and writes it (see read_qpiq/write_qpiq)
  void form_qpiq(const std::string &ident = "");

  //! Forms Q*Pi*Q (with screening, if set) from a given polarisation operator
  //! pi_wk (from polarisation_wk()), which may be shared between Feynman
  //! objects that differ only in screening. Must be on the same sub-grid and
  //! frequency grid
  void form_qpiq(const std::vector<std::vector<ComplexRMatrix>> &pi_wk);

  //! True once Q*Pi*Q has been formed
  bool has_qpiq() const { return m_qpiq_wk.size() != 0; }

  //! Reads Q*Pi*Q from the disk cache for ident (file prefix); false if no
  //! matching file (or ident is "" or "false")
  bool read_qpiq(const std::string &ident);

  //! Writes Q*Pi*Q to the disk cache for ident; false if not written
  bool write_qpiq(const std::string &ident);

  //! Calculate Direct part of correlation potential
  GMatrix Sigma_direct(int kappa_v, double en_v,
                       std::optional<int> k = {}) const;

  //! Exchange part of the correlation potential (second order, no
  //! screening), by frequency integration.
  //! @details The first of the two frequency integrals is done analytically
  //! (closing the contour on the core poles), leaving a single integral
  //! along w = omre + iu, on the same grid as the direct term. Cf.
  //! Goldstone::Sigma_exchange
  GMatrix Sigma_exchange(int kappa_v, double en_v) const;

  //! Direct part of correlation potential, for each multipole k separately
  //! (not summed: Sigma_d = sum of these). Same cost as a single Sigma_direct
  //! call, since the Green's functions are shared by all k
  std::vector<GMatrix> Sigma_direct_each_k(int kappa_v, double en_v) const;

  //! Returns (reference to) q^k (radial) matrix. Note: includes drj? No?
  const ComplexRMatrix &get_qk(int k) const { return m_qk.at(std::size_t(k)); }

  //! Returns (ref to) radial exchange matrix Vx_kappa. Nb: includes dri*drj
  const GMatrix &get_Vx_kappa(int kappa) const {
    return m_Vx_kappa.at(Angular::kappa_to_kindex(kappa));
  }

  //! Returns (ref to) summed core projection operator, P = sum_a |a><a|,
  //! for given kappa. No integration measure.
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

  // Components of Fa on the sub-grid: F[0] = f, F[1] = g (if include_G)
  std::vector<std::vector<double>>
  subgrid_components(const DiracSpinor &Fa) const;

  // The two terms of Gamma (see Sigma_exchange), for one core state a, with
  // the Coulomb line q^l_i2 and the angular factor L attached, summed over
  // the internal partial wave and l. p^a is rank one, so the factor of p^a
  // outside the i sum, F_a^mu(r_1) [pa_gex] or F_a^t(r_j) [gex_pa], is left
  // out. Indexed [spinor index](k, gamma):
  //   pa_gex[t](k, gamma)_j2 = sum_{beta l} L^{kl}_{v beta a gamma} [gex^beta(e_a+w) F_a q^l]^t_j2
  //   gex_pa[mu](k, gamma)_12 = sum_{alpha l} L^{kl}_{v a alpha gamma} [gex^alpha(e_a-w) F_a q^l]^mu_12
  struct GammaQ {
    std::vector<LinAlg::Matrix<ComplexRMatrix>> pa_gex, gex_pa;
  };
  GammaQ exchange_Gamma_q(int kappa_v, const DiracSpinor &Fa,
                          std::complex<double> w,
                          const Angular::SixJTable &sixj) const;

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

//! Angular factor for the exchange correlation potential,
//! L^{kl}_{v beta alpha gamma} of Methods Eq. (RadialSigmaExch), including
//! that equation's (-1)^(k+l) / [j_v]:
//!   L = (-1)^(k+l) C~^k_{v alpha} C~^k_{beta gamma} C~^l_{v gamma}
//!       C~^l_{beta alpha} {j_v j_alpha k; j_beta j_gamma l} / [j_v]
//! @details alpha, beta, gamma are the partial waves of the lines 1i, ij, j2
//! (leaving r_1, joining the two internal points, entering r_2); k, l the
//! multipoles of the Coulomb lines 1j, i2. The symmetric C~
//! (Angular::tildeCk_kk) is used: the product of four C^k of Methods is
//! identical, as the (-1)^(j+1/2) phases cancel in pairs.
double L_exchange(int k, int l, int kappa_v, int kappa_alpha, int kappa_beta,
                  int kappa_gamma, const Angular::SixJTable &sixj);

//! Best real part of the frequency, omre, for the direct diagram: the value
//! furthest from any pole of the u=0 integrand.
//! @details With Delta = e_lowest_excited - e_core_max, the window
//! (-Delta, 0) is free of true poles; inside it lie only the weak fictitious
//! poles at core-core energy differences (imperfect core subtraction in the
//! polarisation-loop Gex) and, for non-lowest valence states, the small
//! valence-valence differences (valence-line G). Only the lowest few
//! physical excited states enter, so the valence list suffices: no
//! basis/spectrum needed (a typical valence energy is assumed if the list is
//! empty). Returns the midpoint of the largest pole-free gap; one omre
//! serves all valence states. Set print=true to show Delta, the in-window
//! poles, and the chosen omre.
double best_omre(const std::vector<DiracSpinor> &core,
                 const std::vector<DiracSpinor> &valence, bool print = false);

} // namespace MBPT