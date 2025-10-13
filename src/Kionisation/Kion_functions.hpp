#pragma once
#include "DiracOperator/Operators/jL.hpp"
#include "LinAlg/Matrix.hpp"
class DiracSpinor;
class Grid;
namespace HF {
class HartreeFock;
}

namespace Kion {

//! Methods for treating continuum part;
/*! @details
 Standard is relativistic Hartree-Fock.
 Zeff is H-like with effective Z.
 Approx is step-function: K(E,q) = Theta[E-enk]*K_nk(q)
*/
enum class Method { HF, RPA0, RPA, Zeff, Approx, Error };

//! DM-electron couplings: Vector and PseudoVector are spin-independent only
enum class Coupling { Vector, Scalar, PseudoVector, PseudoScalar, Error };

//! Format for output file.
/*! @details
xyz: For easy 2D interpolation. list formmated with each row 'E q K(E,q)'
gnuplot: For easy plotting. Each column is new E.
matrix: Outputs entire matrix in table form. E and q grids printed prior.
*/
enum class OutputFormat { gnuplot, matrix, xyz, Error };

//! Units used in output file
/*!
 Atomic: [q] = [1/a_0], [E] = Hartree; \n
 Particle: [q] = MeV, [E] = keV;
*/
enum class Units { Atomic, Particle, Error };

//! Checks if radial grid is dense enough at large r for continuum state
bool check_radial_grid(double Emax, double qmax, const Grid &rgrid);

//! Calculates ionisation factor K(E,q) for given core state, Fnk, using
//! standard method. Stored as matrix. use_rpa0 is flag for including
//! lowest-order RPA (i.e., with zero iterations)
LinAlg::Matrix<double>
calculateK_nk(const HF::HartreeFock *vHF, const DiracSpinor &Fnk, int max_L,
              const Grid &Egrid, const DiracOperator::jL *jl,
              bool force_rescale, bool hole_particle, bool force_orthog,
              bool zeff_cont, bool use_rpa0 = false,
              const std::vector<DiracSpinor> &basis = {},
              double ec_cut = 1000.0);

//! Calculates ionisation factor K(E,q), for all core states, in RPA approximation.
//! Uses all-orders RPA, so is quite slow (RPA must be solved for each L and q)
std::vector<LinAlg::Matrix<double>> calculateK_nk_rpa(
    const HF::HartreeFock *vHF, const std::vector<DiracSpinor> &core, int max_L,
    const Grid &Egrid, DiracOperator::jL *jl, bool force_rescale,
    bool hole_particle, bool force_orthog,
    const std::vector<DiracSpinor> &basis, const std::string &atom);

//! Calculates approxinate ionisation factor K_nk(q) for each core state.
//! Stores in matrix K[nk,q].
/*! @details
Then, K(E,q) = Theta[E-enk]*K_nk(q).
K_nk(q) is calculated with continuum energy just above threshold.
*/
LinAlg::Matrix<double>
calculateK_nk_approx(const HF::HartreeFock *vHF,
                     const std::vector<DiracSpinor> &core, int max_L,
                     const DiracOperator::jL *jl, bool force_rescale,
                     bool hole_particle, bool force_orthog, bool zeff_cont,
                     bool use_rpa, const std::vector<DiracSpinor> &basis);

//! Converts K_nk(q) to K(E,q) - for ease of plotting/comparison
LinAlg::Matrix<double>
convert_K_nk_approx_to_std(const LinAlg::Matrix<double> &Kaprx,
                           const Grid &Egrid,
                           const std::vector<DiracSpinor> &core);

//! Writes ouput file in matrix form
/*! @details
matrix : Outputs entire matrix in table form. E and q grids printed prior.
In K[E,q] form: each column is different q
*/
void write_to_file_matrix(const LinAlg::Matrix<double> &K, const Grid &E_grid,
                          const Grid &q_grid, const std::string &filename,
                          int num_digits = 5, Units units = Units::Particle);

//! Writes ouput file in 'xyz' form: for easy 2D interpolation
/*! @details
xyz: For easy 2D interpolation. list formmated with each row 'E q K(E,q)'
*/
void write_to_file_xyz(const LinAlg::Matrix<double> &K, const Grid &E_grid,
                       const Grid &q_grid, const std::string &filename,
                       int num_digits = 5, Units units = Units::Particle);

//! Writes ouput file in 'gnuplot' form: for easy plotting
/*! @details
gnuplot: For easy plotting. Each column is new E.
First column is q values.
First row is E values (use will `t columnheader(i)`)
*/
void write_to_file_gnuplot(const LinAlg::Matrix<double> &K, const Grid &E_grid,
                           const Grid &q_grid, const std::string &filename,
                           int num_digits = 5, Units units = Units::Particle);

//! Writes to file: formats is a vector of formats, will write to each listed format
void write_to_file(const std::vector<OutputFormat> &formats,
                   const LinAlg::Matrix<double> &K, const Grid &E_grid,
                   const Grid &q_grid, const std::string &filename,
                   int num_digits = 5, Units units = Units::Particle);

//! Writes the 'approximate tables' to file. Each column is a new state, K_nk(q)
void write_approxTable_to_file(const LinAlg::Matrix<double> &K,
                               const std::vector<DiracSpinor> &core,
                               const Grid &q_grid, const std::string &filename,
                               int num_digits, Units units);

} // namespace Kion