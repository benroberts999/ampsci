#pragma once
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/Ladder.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Physics/AtomData.hpp" //DiracSEnken
#include "Wavefunction/DiracSpinor.hpp"
#include <cassert>
#include <vector>
class Grid;
namespace HF {
class HartreeFock;
}

//! Many-body perturbation theory
namespace MBPT {

//******************************************************************************
// Helper classes (store parameters):

enum class Method { Goldstone, Feynman };

struct Sigma_params {
  Method method;
  int min_n_core;
  bool include_G{false};

  // Following only for Feynman method
  int max_l_excited{6};
  bool GreenBasis{false};
  bool PolBasis{false};
  double real_omega{-0.2};
  double w0{0.01};
  double w_ratio{1.5};
  bool screenCoulomb{false};
  bool holeParticle{false};
  std::string ladder_file{""};

  std::vector<double> fk{};   // this for Goldstone too..
  std::vector<double> etak{}; // this for Goldstone too..
};

struct rgrid_params {
  double r0;
  double rmax;
  std::size_t stride;
};

struct wgrid_params {
  // Only used in Feynman method:
  double w0;
  double wmax;
  double ratio;
  // XX add exchange stride ?
};

//******************************************************************************
/*!
@brief Pure virtual. Calculates + stores Correlation Potential operator
@details
Takes core+excited basis orbitals, and a list of energies.
Must be one energy per kappa, in order starting from k=-1 (no gaps).
Correlation potential is formed for each kappa for which an energy is given, at
the given energy.

Can calculate second-order MBPT energy shifts without forming matrix (give blank
energy list).

in_fname: will read/write to this file.
If the file exists, will read Sigma from it. If file doesn't exist, will
calculate Sigma and write it to file. If blank filename is given, will not try
to write.
 - Note: program appends ".Sigma" extension, don't include this in filename.
 - Note: only writes sigma, not scaling factors

*/
class CorrelationPotential {
protected:
  //! If given (optional) en_list, will form Sigma matrix.
  //! @details if given filename, will read Sigma in from file instead of
  //! calculating it. If given filename for file that doesn't exist yet,
  //! will write sigma to that file.
  CorrelationPotential(const HF::HartreeFock *const in_hf,
                       const std::vector<DiracSpinor> &basis,
                       const Sigma_params &sigp, const rgrid_params &subgridp);

  // Try to "undelete" these??
  CorrelationPotential &operator=(const CorrelationPotential &) = delete;

public:
  CorrelationPotential(const CorrelationPotential &) = default;

public:
  virtual ~CorrelationPotential() = default;

  // Calculates Sigma, for given kappa, energy, and stores.
  // Also stored "lookup table" including n
  // n may be zero (means)
  // This should be virtual?
  // nb: should do nothing if sigma already exists??
  virtual void formSigma(int kappa, double en, int n = 0) {
    (void)kappa;
    (void)en;
    (void)n; // don't warn on unsused, want named
    assert(false && "Cannot call formSigma on copied CorrelationPotential!");
  };

  const GMatrix *getSigma(int n, int kappa) const;

  //! returns Spinor: Sigma|Fv>
  //! @details If Sigma for kappa_v doesn't exist, returns |0>. m_Sigma_kappa
  //! calculated at the energy given in 'form_Sigma' (or on construct)
  DiracSpinor SigmaFv(const DiracSpinor &Fv) const;
  DiracSpinor operator()(const DiracSpinor &Fv) const { return SigmaFv(Fv); }

protected:
  // // n=0 means get Sigma for lowest available n
  std::size_t getSigmaIndex(int n, int kappa) const;

  // //! Forms Correlation potential matrix [called on construct]
  // //! @details Forms for each kappa at given energy. One energy must be given
  // //! per kappa.
  // void form_Sigma(const std::vector<DiracSpinor> &valence = {},
  //                 const std::string &out_fname = "");

public:
  //! Stores scaling factors, lambda, for each kappa (Sigma -> lamda*Sigma)
  void scale_Sigma(const std::vector<double> &lambda_kappa) {
    m_lambda_kappa = lambda_kappa;
  }

  void scale_Sigma(int n, int kappa, double lambda);

  //! Prints the scaling factors to screen
  void print_scaling() const;
  //! Prints the sub-grid parameters to screen
  void print_subGrid() const;

  void print_info() const {
    print_subGrid();
    if (!m_nk.empty())
      std::cout << "Have Sigma for:\n";
    for (const auto [n, k, en] : m_nk) {
      std::cout << n << " " << AtomData::kappa_symbol(k) << " en=" << en
                << "\n";
    }
  }

  bool empty() const { return m_Sigma_kappa.empty(); }

  //! Calculates <Fv|Sigma|Fw> from scratch, at Fv energy [full grid + fg+gg]
  //! @details Note: uses basis, so if reading Sigma from file, and no basis
  //! given, will return all 0.0
  double SOEnergyShift(const DiracSpinor &Fv, const DiracSpinor &Fw,
                       int max_l = 99) const;

  static DiracSpinor Sigmal_Fv(const DiracSpinor &v, const Coulomb::YkTable &yk,
                               const Coulomb::LkTable &lk,
                               const std::vector<DiracSpinor> &core,
                               const std::vector<DiracSpinor> &excited);

  RDMatrix<double> Sigma_l(const DiracSpinor &v, const Coulomb::YkTable &qk,
                           const Coulomb::LkTable &lk,
                           const std::vector<DiracSpinor> &core,
                           const std::vector<DiracSpinor> &excited) const;

  int maxk() const { return m_maxk; }

  // make virtual??
  // Read and writes Sigma (G) matrix to file
  bool read_write(const std::string &fname, IO::FRW::RoW rw);

  // Acts Gmatirx (G) matrix onto Fv. Interpolates from sub-grid
  DiracSpinor act_G_Fv(const GMatrix &Gmat, const DiracSpinor &Fv) const;
  double act_G_Fv_2(const DiracSpinor &Fa, const GMatrix &Gmat,
                    const DiracSpinor &Fb) const;

protected:
  void setup_subGrid(double rmin, double rmax);

  // Adds new |ket><bra| term to G; uses sub-grid
  void addto_G(GMatrix *Gmat, const DiracSpinor &ket, const DiracSpinor &bra,
               const double f = 1.0) const;

  std::vector<DiracSpinor> copy_holes(const std::vector<DiracSpinor> &basis,
                                      const std::vector<DiracSpinor> &core,
                                      int n_min) const;
  std::vector<DiracSpinor>
  copy_excited(const std::vector<DiracSpinor> &basis,
               const std::vector<DiracSpinor> &core) const;

public:
  // Project subgrid onto full grid
  // std::size_t ri_subToFull(std::size_t i) const;
  // get value of dr on subgrid
  // double dr_subToFull(std::size_t i) const;

  std::shared_ptr<const Grid> p_gr;
  // occupied (holes) and excited (virtual) states. Holes includes from n>=nmin
  const std::vector<DiracSpinor> m_holes, m_excited;

protected:
  // Coulumb Y^k_eh (excited, holes) table (includes C^k)
  Coulomb::YkTable m_yeh;
  // maximum multipolarity, k, = 2*max(j) {max 2*j in core/basis}
  int m_maxk;
  // 6j lookup table
  // Angular::SixJ m_6j;
  Angular::SixJTable m_6j;

  const bool m_ratio_ladder_method = false;
  std::string m_ladder_file;
  std::optional<Coulomb::LkTable> m_lk{std::nullopt};

  // SubGrid values:
  // Stride: only include every 'stride' grid points (between min/max)
  std::size_t m_stride;
  // Number of points in subgrid
  std::size_t m_subgrid_points{};
  // First real-grid point included r[m_imin]=r_min {real grid}
  std::size_t m_imin{};
  // Values for r on the subgrid
  std::vector<double> m_subgrid_r{};

  // m_Sigma_kappa: holds Sigma matrix for each partial-wave, kappa
  std::vector<GMatrix> m_Sigma_kappa{};
  // Lambda (fitting factors) for each kappa
  std::vector<double> m_lambda_kappa{};
  std::vector<AtomData::DiracSEnken> m_nk{};

  // Options for sub-grid, and which matrices to include
  const bool m_include_G;

  bool same_as_fileQ{false};

  // Effective screening parameters
  std::vector<double> m_fk{};  // e.g., {0.72, 0.62, 0.83, 0.89, 0.94, 1.0};
  std::vector<double> m_eta{}; // e.g., {1.1, 1.4, 1.0};

  // XXX KILL ?
  double get_fk(int k) const {
    if (k < int(m_fk.size())) {
      return m_fk[std::size_t(k)];
    }
    return 1.0;
  }
  double get_eta(int k) const {
    if (k < int(m_eta.size())) {
      return m_eta[std::size_t(k)];
    }
    return 1.0;
  }
};

} // namespace MBPT
