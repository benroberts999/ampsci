#pragma once
#include "DiracOperator/Operators/jL.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/SphericalBessel.hpp"
#include <array>
#include <vector>
class DiracSpinor;
class Grid;
class Wavefunction;
namespace HF {
class HartreeFock;
}

namespace Kion {

//! DM-electron couplings
enum class Coupling { Vector, Scalar, AxialVector, PseudoScalar, Error };

//! Format for output file.
/*! @details
xyz: For easy 2D interpolation. list formmated with each row 'E q K(E,q)'
gnuplot: For easy plotting. Each column is new E.
matrix: Outputs entire matrix in table form. E and q grids printed prior.
*/
enum class OutputFormat { matrix, xyz, Error };

//! Units used in output file
/*! @details
 Atomic: [q] = [1/a_0], [E] = Hartree;
 Particle: [q] = eV, [E] = eV;
*/
enum class Units { Atomic, Particle, Error };

//! Checks if radial grid is dense enough at large r for continuum state
bool check_radial_grid(double Emax, double qmax, const Grid &rgrid);

//! Calculates all 13 form factors (V,A,S,P) for a single core state Fa.
/*! @details
 Returns an array of 13 matrices: {K_VT, K_VE, K_VM, K_VL, K_T5, K_E5,
 K_M5, K_L5, K_X, K_X5, K_Z, K_S, K_S5}. 
 @note: Matrix will be empty (0x0) if not calculated (set by bool).
 @note: order is important
*/
std::array<LinAlg::Matrix<double>, 13> calculate_formFactors_nk(
  const Wavefunction &wf, const DiracSpinor &Fa, int lc_min, int lc_max,
  int Kmin, int Kmax, const std::vector<double> &Egrid,
  const std::vector<double> &qgrid, bool diagonal_Eq, bool force_rescale,
  bool hole_particle, bool force_orthog, bool low_q,
  const SphericalBessel::JL_table &jK_tab, bool vectorQ, bool spatialQ,
  bool axialQ, bool XvvQ, bool YaaQ, bool ZvaQ, bool scalarQ,
  bool pseudoscalarQ);

//! Calculates ionisation factor K(E,q) for given core state, Fnk, using
//! standard method. Stored as matrix. use_rpa0 is flag for including
//! lowest-order RPA (i.e., with zero iterations)
LinAlg::Matrix<double>
calculateK_nk(const HF::HartreeFock *vHF, const DiracSpinor &Fnk, int max_L,
              const Grid &Egrid, const DiracOperator::jL *jl,
              bool force_rescale, bool hole_particle, bool force_orthog,
              bool zeff_cont, bool zeff_bound, double ec_cut = 1.0e99);

//! Writes ouput file in matrix form
/*! @details
matrix : Outputs entire matrix in table form. E and q grids printed prior.
In K[E,q] form: each column is different q
*/
void write_to_file_matrix(const LinAlg::Matrix<double> &K,
                          const std::vector<double> &E_grid,
                          const std::vector<double> &q_grid,
                          const std::string &filename, int num_digits = 5,
                          Units units = Units::Particle);

//! Writes ouput file in 'xyz' form: for easy 2D interpolation
/*! @details
xyz: For easy 2D interpolation. list formmated with each row 'E q K(E,q)'
*/
void write_to_file_xyz(const std::string &filename,
                       const std::vector<double> &E_grid,
                       const std::vector<double> &q_grid,
                       const std::vector<std::string> &titles,
                       const std::vector<std::string> &descriptions,
                       std::vector<LinAlg::Matrix_view<const double>> factors,
                       Units units = Units::Particle, int num_digits = 6,
                       bool diagonal = false);

void write_to_file_xyz_13(
  const std::string &filename, const std::vector<double> &E_grid,
  const std::vector<double> &q_grid, const std::vector<std::string> &titles,
  const std::vector<std::string> &descriptions,
  const std::array<LinAlg::Matrix<double>, 13> K_factors,
  Units units = Units::Particle, int num_digits = 6, bool diagonal = false);

} // namespace Kion