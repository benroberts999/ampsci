#pragma once
#include "Angular/Angular_tables.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <vector>
class Grid;

//! Many-body perturbation theory
namespace MBPT {

//! Holds Green's fn operator of form: |ket><bra| [4x4 matrix of NxN matrix]
struct GMatrix {
  GMatrix(int in_size);
  int size; // careful, should be const, but then need copy construct?
  LinAlg::SqMatrix ff, fg, gf, gg;
  //! Sets all matrix elements to zero
  void zero();
  //! Can add/subtract matrices (in place)
  GMatrix &operator+=(const GMatrix &rhs);
  GMatrix &operator-=(const GMatrix &rhs);
};

//******************************************************************************
/*!
@brief Calculates + stores Correlation Potential operator (2nd order)
@details
Takes core+excited basis orbitals, and a list of energies.
Should be one energy per kappa, in order starting from k=-1 (no gaps).
Correlation potential is formed for each kappa for which an energy is given, at
the given energy.
Can calculate energy shifts without forming matrix (give blank energy list).

in_fname: will read/write to this file.
If the file exists, will read Sigma from it. If file doesn't exist, will
calculate Sigma and write it to file. If blank filename is given, will not try
to write.
 - Note: program appends ".Sigma" extension, don't include this in filename.
 - Note: only writes sigma, not scaling factors

Sigma stored in a matrix consisting of 4 'G' terms:

Four second-order diagrams:
  - Diagram (a):
\f[ G_a = \frac{|Q^k_amn><Q^k_amn|}{ [k][j] \, de_{amn} } \f]
  - Diagram (b) (exchange):
\f[ G_b = \frac{|Q^k_amn><P^k_amn|}{ [k][j] \, de_{amn} } \f]
  - Diagram (c):
\f[ G_c = \frac{|Q^k_nba><Q^k_nba|}{ [k][j] \, de_{nba} } \f]
  - Diagram (d) (exchange):
\f[ G_d = \frac{|Q^k_nba><P^k_nba|}{ [k][j] \, de_{nba} } \f]
where:
All indeces are summed over,
a & b are core states, n & m are virtual excited states,
k is multipolarity [Coloulmb expansion], and
\f$ de_{xyz} = e_v + e_x - e_y - e_z \f$

*/
class CorrelationPotential {
public:
  //! If given (optional) en_list, will form Sigma matrix.
  //! @details if given filename, will read Sigma in from file instead of
  //! calculating it. If given filename for file that doesn't exist yet, will
  //! write sigma to that file.
  CorrelationPotential(const Grid &gr,
                       const std::vector<DiracSpinor> &core = {},
                       const std::vector<DiracSpinor> &excited = {},
                       const int in_stride = 4,
                       const std::vector<double> &en_list = {},
                       const std::string &in_fname = "");

  CorrelationPotential &operator=(const CorrelationPotential &) = delete;
  CorrelationPotential(const CorrelationPotential &) = delete;
  ~CorrelationPotential() = default;

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

  //! returns Spinor: Sigma|Fv>
  //! @details If Sigma for kappa_v doesn't exist, returns |0>. Sigma_kappa
  //! calculated at the energy given in 'form_Sigma' (or on construct)
  DiracSpinor operator()(const DiracSpinor &Fv) const { return Sigma2Fv(Fv); }

  //! Calculates <Fv|Sigma|Fw> from scratch, at Fv energy [full grid + fg+gg]
  //! @details Note: uses basis, so if reading Sigma from file, and no basis
  //! given, will return all 0.0
  double operator()(const DiracSpinor &Fv, const DiracSpinor &Fw) const {
    return Sigma2vw(Fv, Fw);
  }

private:
  DiracSpinor Sigma2Fv(const DiracSpinor &Fv) const;
  double Sigma2vw(const DiracSpinor &Fv, const DiracSpinor &Fw) const;
  void setup_subGrid();

private:
  // main routine, filles G (Sigma) matrix using basis
  void fill_Gkappa(GMatrix *Gmat, const int kappa, const double en);
  // Adds new |ket><bra| term to G; uses sub-grid
  void addto_G(GMatrix *Gmat, const DiracSpinor &ket, const DiracSpinor &bra,
               const double f = 1.0) const;
  // Acts Sigma (G) matrix onto Fv. Interpolates from sub-grid
  DiracSpinor Sigma_G_Fv(const GMatrix &Gmat, const DiracSpinor &Fv) const;
  // Read and writes Sigma (G) matrix to file
  void read_write(const std::string &fname, IO::FRW::RoW rw);

private:
  const Grid *const p_gr;
  const std::vector<DiracSpinor> m_core;
  const std::vector<DiracSpinor> m_excited;
  Coulomb::YkTable m_yec; // constains Ck and Y_ec(r)
  const int m_maxk;
  Angular::SixJ m_6j;
  int stride;

  int stride_points{};
  int imin{};
  std::vector<double> r_stride{};
  std::vector<GMatrix> G_kappa{};

  std::vector<double> m_lambda_kappa{};

  // Options for sub-grid, and which matrices to include
  static constexpr bool include_FG = false;
  static constexpr bool include_GG = false;
};

} // namespace MBPT