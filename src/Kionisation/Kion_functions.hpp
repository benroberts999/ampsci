#pragma once
#include "DiracOperator/Operators/jL.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <array>
#include <cmath>
#include <string>
#include <vector>
class DiracSpinor;
class Grid;
class Wavefunction;
namespace HF {
class HartreeFock;
}

//! Functions for atomic ionisation form factors
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

/*! @brief
  Method used to solve bound/continuum states for form factors

  @details
  HF: real (Hartree-Fock) bound and continuum states (standard method).
  Zeff: H-like (Zeff) bound and continuum states, solved numerically
  with DiracODE.
  ZeffAnalytic: H-like (Zeff) bound and continuum states, using exact
  analytic Dirac-Coulomb functions. Relativistic continuum requires FLINT
  (see DiracContinuum::available).
*/
enum class AtomicMethod { HF, Zeff, ZeffAnalytic };

//! Parses string (HF, Zeff, ZeffAnalytic) to AtomicMethod (case-insensitive).
//! Unknown input: warns, defaults to HF.
AtomicMethod parseStatesMethod(const std::string &in_method);
//! StatesMethod to string (HF, Zeff, ZeffAnalytic)
std::string parseStatesMethod(const AtomicMethod &in_method);

//! Effective charge from binding energy: Zeff = n * sqrt(-2*en).
//! Same Zeff as used by DarkARC (see arxiv:1912.08204).
inline double Zeff_real(double en, int n) {
  return n * std::sqrt(std::abs(2.0 * en));
}

//! Checks if radial grid is dense enough at large r for continuum state,
//! and (roughly) for the maximum safe q.
bool check_radial_grid(double Emax, double qmax, const Grid &rgrid,
                       double alpha = PhysConst::alpha);

/*! 
  @brief Calculates all 13 form factors (V,A,S,P) for a single core state Fa.

  @details
  Returns an array of 13 matrices: {K_VT, K_VE, K_VM, K_VL, K_T5, K_E5,
  K_M5, K_L5, K_X, K_X5, K_Z, K_S, K_S5}.
  @note: Matrix will be empty (0x0) if not calculated (set by bool).
  @note: order is important

  Optionally (method != AtomicMethod::HF), uses H-like (Zeff) states for
  both the bound state and the continuum, solved either numerically
  (DiracODE) or with exact analytic Dirac-Coulomb functions (see
  AtomicMethod). Zeff is zeff_constant if non-zero, else the "real"
  Zeff = n*sqrt(-2*en) from the binding energy (see Zeff_real). Continuum
  energies (ec = E + en) and occupation always use the real (input) Fa.
  force_rescale and hole_particle have no effect for Zeff states.
*/
std::array<LinAlg::Matrix<double>, 13> calculate_formFactors_nk(
  const HF::HartreeFock *vHF, const DiracSpinor &Fa, int lc_min, int lc_max,
  double ec_min, double ec_max, bool force_rescale, bool hole_particle,
  bool force_orthog, const std::vector<double> &Egrid,
  const std::vector<double> &qgrid, bool diagonal_Eq, bool low_q, bool high_q,
  const SphericalBessel::JL_table &jK_tab, int Kmin, int Kmax, bool vectorQ,
  bool axialQ, bool scalarQ, bool pseudoscalarQ, bool spatialQ,
  AtomicMethod method = AtomicMethod::HF, double zeff_constant = 0.0);

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