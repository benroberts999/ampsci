#pragma once
#include "Angular/Angular_tables.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/GreenMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>
class Grid;
namespace HF {
class HartreeFock;
}

//! Many-body perturbation theory
namespace MBPT {

//******************************************************************************
// Helper function:
inline int find_max_tj(const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &excited) {
  // returns maximum value of 2*j in {core,excited}
  auto maxtj1 = core.empty() ? 0
                             : std::max_element(core.cbegin(), core.cend(),
                                                DiracSpinor::comp_j)
                                   ->twoj();
  auto maxtj2 = excited.empty()
                    ? 0
                    : std::max_element(excited.cbegin(), excited.cend(),
                                       DiracSpinor::comp_j)
                          ->twoj();
  return std::max(maxtj1, maxtj2);
}
inline int find_max_l(const std::vector<DiracSpinor> &orbs) {
  return orbs.empty()
             ? 0
             : std::max_element(orbs.cbegin(), orbs.cend(), DiracSpinor::comp_l)
                   ->l();
}
//******************************************************************************

// using GMatrix = GreenMatrix<LinAlg::SqMatrix>;
// using ComplexGMatrix = GreenMatrix<LinAlg::ComplexSqMatrix>;
// using ComplexDouble = LinAlg::ComplexDouble;

enum class Method { Goldstone, Feynman };

struct Sigma_params {
  Method method;
  int min_n_core;
  // Following only for Feynman method
  int max_l_excited;
  bool GreenBasis;
  bool PolBasis;
  double real_omega;
  bool screenCoulomb;
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

  CorrelationPotential &operator=(const CorrelationPotential &) = delete;
  CorrelationPotential(const CorrelationPotential &) = delete;

public:
  virtual ~CorrelationPotential() = default;

  //! Forms Correlation potential matrix [called on construct]
  //! @details Forms for each kappa at given energy. One energy must be given
  //! per kappa.
  void form_Sigma(const std::vector<double> &en_list = {},
                  const std::string &out_fname = "");

  //! Stores scaling factors, lambda, for each kappa (Sigma -> lamda*Sigma)
  void scale_Sigma(const std::vector<double> &lambda_kappa) {
    m_lambda_kappa = lambda_kappa;
  }

  //! Prints the scaling factors to screen
  void print_scaling() const;
  //! Prints the sub-grid parameters to screen
  void print_subGrid() const;

  //! returns Spinor: Sigma|Fv>
  //! @details If Sigma for kappa_v doesn't exist, returns |0>. m_Sigma_kappa
  //! calculated at the energy given in 'form_Sigma' (or on construct)
  DiracSpinor SigmaFv(const DiracSpinor &Fv) const;
  DiracSpinor operator()(const DiracSpinor &Fv) const { return SigmaFv(Fv); }

  //! Calculates <Fv|Sigma|Fw> from scratch, at Fv energy [full grid + fg+gg]
  //! @details Note: uses basis, so if reading Sigma from file, and no basis
  //! given, will return all 0.0
  double SOEnergyShift(const DiracSpinor &Fv, const DiracSpinor &Fw,
                       int max_l = 99) const;

public:
  void setup_subGrid(double rmin, double rmax);

  // main routine, filles Sigma matrix using basis [Goldstone]
  virtual void fill_Sigma_k(GMatrix *Gmat, const int kappa,
                            const double en) = 0;
  // Fills Sigma matrix using Feynman technique
  // void fill_Sigma_k_Feyn(GMatrix *Gmat, const int kappa, const double en);
  // Adds new |ket><bra| term to G; uses sub-grid
  void addto_G(GMatrix *Gmat, const DiracSpinor &ket, const DiracSpinor &bra,
               const double f = 1.0) const;

  // Acts Sigma (G) matrix onto Fv. Interpolates from sub-grid
  DiracSpinor Sigma_G_Fv(const GMatrix &Gmat, const DiracSpinor &Fv) const;

  // Read and writes Sigma (G) matrix to file
  bool read_write(const std::string &fname, IO::FRW::RoW rw);

  // Project subgrid onto full grid
  std::size_t ri_subToFull(std::size_t i) const;
  // get value of dr on subgrid
  double dr_subToFull(std::size_t i) const;

  std::vector<DiracSpinor> copy_holes(const std::vector<DiracSpinor> &basis,
                                      const std::vector<DiracSpinor> &core,
                                      int n_min) const;
  std::vector<DiracSpinor>
  copy_excited(const std::vector<DiracSpinor> &basis,
               const std::vector<DiracSpinor> &core) const;

protected:
  std::shared_ptr<const Grid> p_gr;
  // occupied (holes) and excited (virtual) states. Holes includes from n>=nmin
  const std::vector<DiracSpinor> m_holes, m_excited;
  // Coulumb Y^k_eh (excited, holes) table (includes C^k)
  Coulomb::YkTable m_yeh;
  // maximum multipolarity, k, = 2*max(j) {max 2*j in core/basis}
  int m_maxk;
  // 6j lookup table
  Angular::SixJ m_6j;

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

  // Options for sub-grid, and which matrices to include
  static constexpr bool m_include_G = false;
};

} // namespace MBPT
