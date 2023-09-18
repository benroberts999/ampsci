#pragma once
#include "CorrelationPotential.hpp" // for rgrid_params
#include "HF/HartreeFock.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <vector>

namespace MBPT {

//! Options for including hole-particle interaction. include mean all k;
//! include_k0 means k=0 term only
enum class HoleParticle { exclude, include, include_k0 };

//! Options for including Screening
enum class Screening { exclude, include, only };

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
  std::size_t m_i0, m_imax, m_stride, m_subgrid_points;

  // maximum kappa_index appearing in the core
  int m_max_ki_core;
  // maximum kappa index to include for internal lines
  int m_max_ki;
  // maximum multipolarity, k
  int m_max_k;
  // Lowest n to polarise in polarisation operator
  int m_min_core_n;
  // real part of frequency for integration
  double m_omre;
  // Stores imaginary frequency grid (setup post-construction)
  Grid m_wgrid;

  // Option to print to screen, or calculate quietly
  bool m_verbose;

  // Option: include hole-particle interaction into polarisation operator
  bool m_hole_particle;
  // Option to include higher-k in hp interaction
  bool m_include_higher_order_hp;
  // Include screening correction to Coulomb line in Q*Pi*Q => Q*Pi*X*Q
  bool m_screen_Coulomb;
  // _only_ include screening correction in Q*Pi*Q => Q*Pi*[X-1]*Q
  bool m_only_screen;

  // Coulomb operator; include right integration measure
  // q_ij := (r>^k/r<^{k+1}) * dr_j
  std::vector<ComplexGMatrix> m_qk{};
  // Left integration measure:
  ComplexGMatrix m_dri;
  // Right integration measure:
  ComplexGMatrix m_drj;

  // Core projection operators (one for each core state)
  std::vector<ComplexGMatrix> m_pa{};
  // Hartree-Fock Exchange matrix (one for each happa)
  std::vector<GMatrix> m_Vx_kappa{};

  // Effective Q*Pi*Q operator: for each imaginary omega, and each k
  std::vector<std::vector<ComplexGMatrix>> m_qpiq_wk{};

  // For now, just for testing: switch between complex green methods
  bool m_Complex_green_method = false;

public:
  //! Construct Feynman diagram
  //! @details rgrid_params={r0, rmax, stride}; omre is real part of frequency,
  //! w_0 is initial point along imaginary freq. axis;
  //! w_ratio is ratio used for logarithic omega grid (integration);
  //! scr_option and hp_option are screening and hole-particle interactions;
  //! max_l is maximum l to include for internal lines (Green's functions);
  //! n_min_core is minimum n to include in polarisation loop **
  //! ** Currently have issue: polarising deep n leads to failure?
  Feynman(const HF::HartreeFock *vHF, MBPT::rgrid_params subgrid, double omre,
          double w_0, double w_ratio, Screening scr_option,
          HoleParticle hp_option, int max_l = 6, int n_min_core = 1,
          bool verbose = true);

public:
  //! Calculates Green's function for kappa, and complex energy
  ComplexGMatrix green(int kappa, std::complex<double> en,
                       GreenStates states) const;

  //! Polarisation operator pi^k(w), for multipolarity k
  ComplexGMatrix polarisation_k(int k, std::complex<double> omega) const;

  //! Calculate Direct part of correlation potential
  GMatrix Sigma_direct(int kappa_v, double en_v,
                       std::optional<int> k = {}) const;

  //! Returns (reference to) q^k (radial) matrix. Note: includes drj
  const ComplexGMatrix &get_qk(int k) const { return m_qk.at(std::size_t(k)); }

  //! Returns (ref to) radial exchange matrix Vx_kappa. Nb: includes dri*drj
  const GMatrix &get_Vx_kappa(int kappa) const {
    return m_Vx_kappa.at(std::size_t(Angular::indexFromKappa(kappa)));
  }

  //! Returns stride used for sub-grid
  std::size_t stride() const { return m_stride; }

private:
  // hack: checks if n_min_core is OK, updates if not
  void check_min_n();
  // forms Qk matrices, as well as dri, drj
  void form_qk();
  // Forms core projection operators
  void form_pa();
  // Forms HF exchange potential matrix
  void form_vx();
  // Sets up imaginary frequency grid
  Grid form_w_grid(double w0, double wratio) const;
  // Constructs the Q*Pi*Q Matrix along w grid, for each k
  void form_qpiq();

  // Screening factor X = [1 + i qk*pik]^-1
  ComplexGMatrix X_screen(const ComplexGMatrix &pik,
                          const ComplexGMatrix &qk) const;

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
  ComplexGMatrix green_hf_v2(int kappa, std::complex<double> en,
                             const DiracSpinor *Fc_hp = nullptr) const;

  // Given Gr(wr), and wi (wr is real), returns Gr(wr + i*wi)
  ComplexGMatrix green_to_complex(const ComplexGMatrix &Gr,
                                  double om_imag) const;

  // Green's function, only due to core states
  ComplexGMatrix green_core(int kappa, std::complex<double> en) const;

  // Greens function, only due to excited states (by orthog to core)
  ComplexGMatrix green_excited(int kappa, std::complex<double> en,
                               const DiracSpinor *Fc_hp = nullptr) const;

  // Forces Green's fn to be orthog. to core
  ComplexGMatrix orthogonalise_wrt_core(const ComplexGMatrix &g_in,
                                        int kappa) const;

  // Given Dirac solutions regular at 0 (x0) and infinity (xI), forms "local"
  // Green's function
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

  [[nodiscard]] const ComplexGMatrix &get_dri() const { return m_dri; }
  [[nodiscard]] const ComplexGMatrix &get_drj() const { return m_drj; }

public:
  Feynman &operator=(const Feynman &) = default;
  Feynman(const Feynman &) = default;
  ~Feynman() = default;
};

} // namespace MBPT