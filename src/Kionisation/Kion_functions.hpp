#pragma once
#include "DiracOperator/Operators/jL.hpp"
#include "LinAlg/Matrix.hpp"
class DiracSpinor;
class Grid;
namespace HF {
class HartreeFock;
}

namespace Kion {

//! DM-electron couplings: Vector and PseudoVector are spin-independent only
enum class Coupling { Vector, Scalar, PseudoVector, PseudoScalar, Error };

//! Format for output file.
/*! @details
xyz: For easy 2D interpolation. list formmated with each row 'E q K(E,q)'
gnuplot: For easy plotting. Each column is new E.
matrix: Outputs entire matrix in table form. E and q grids printed prior.
*/
enum class OutputFormat { gnuplot, gnuplot_E, matrix, xyz, Error };

//! Units used in output file
/*! @details
 Atomic: [q] = [1/a_0], [E] = Hartree;
 Particle: [q] = eV, [E] = eV;
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
              bool zeff_cont, double ec_cut = 1.0e99);

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
void write_to_file_xyz(const LinAlg::Matrix<double> &K,
                       const std::vector<double> &E_grid,
                       const std::vector<double> &q_grid,
                       const std::string &filename, int num_digits = 5,
                       Units units = Units::Particle);

//! Order is important!
void write_to_file_xyz_set(const std::vector<double> &E_grid,
                           const std::vector<double> &q_grid,
                           const std::string &filename, int num_digits,
                           Units units, const LinAlg::Matrix<double> &K_T = {},
                           const LinAlg::Matrix<double> &K_E = {},
                           const LinAlg::Matrix<double> &K_M = {},
                           const LinAlg::Matrix<double> &K_L = {},
                           const LinAlg::Matrix<double> &K_X = {});

//! Writes ouput file in 'gnuplot' form: for easy plotting as function of q
/*! @details
gnuplot: For easy plotting. Each column is new E.
First column is q values.
First row is E values (use with `t columnheader(i)`)
*/
void write_to_file_gnuplot(const LinAlg::Matrix<double> &K,
                           const std::vector<double> &E_grid,
                           const std::vector<double> &q_grid,
                           const std::string &filename, int num_digits = 5,
                           Units units = Units::Particle);

//! Writes ouput file in 'gnuplot' form: for easy plotting as function of E.
/*! @details
gnuplot: For easy plotting. Each column is new q.
First column is E values.
First row is q values (use with `t columnheader(i)`)
*/
void write_to_file_gnuplot_E(const LinAlg::Matrix<double> &K,
                             const std::vector<double> &E_grid,
                             const std::vector<double> &q_grid,
                             const std::string &filename, int num_digits = 5,
                             Units units = Units::Particle);

//! Writes to file: formats is a vector of formats, will write to each listed format
void write_to_file(const std::vector<OutputFormat> &formats,
                   const LinAlg::Matrix<double> &K,
                   const std::vector<double> &E_grid,
                   const std::vector<double> &q_grid,
                   const std::string &filename, int num_digits = 5,
                   Units units = Units::Particle);

} // namespace Kion