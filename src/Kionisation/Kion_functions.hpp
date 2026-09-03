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
  const std::vector<double> &qgrid, bool diagonal_Eq, bool low_q,
  const SphericalBessel::JL_table &jK_tab, int Kmin, int Kmax, bool vectorQ,
  bool axialQ, bool scalarQ, bool pseudoscalarQ, bool spatialQ,
  AtomicMethod method = AtomicMethod::HF, double zeff_constant = 0.0);

/*!
  @brief Options for the core-polarisation (continuum RPA) dressing of the
  form factors; see @ref calculate_formFactors_rpa.
*/
struct RPAOptions {
  //! Maximum TDHF iterations per driven solve, and per K-matrix column
  int max_its{60};
  //! Convergence target of the TDHF iterations
  double eps{1.0e-5};
  //! Unitarise: physical (outgoing-wave) amplitudes via the on-shell K
  //! matrix, at the cost of one extra field-free solve per open channel per
  //! (E, rank, parity). If false, the standing-wave amplitudes are used
  //! directly: these have real poles just above thresholds, and need the
  //! HF continuum bra of each orbital.
  bool unitarise{true};
  //! Energy deposits above this (atomic units) are not dressed: their rows
  //! are left zero. Core polarisation is a low-energy effect, while the
  //! solve is at its most expensive when every shell is open. Zero or
  //! negative: no limit.
  double E_max{0.0};
  //! Parallelise the driven solves over momentum transfer q (one TDHF
  //! instance, with its own Anderson history, per thread) rather than over
  //! the channel tasks inside each solve. Better load balance; memory
  //! scales with the thread count (see ExternalField::TDHFcntm).
  bool parallel_q{true};
};

/*!
  @brief Form factors with core polarisation (RPA) for every core orbital,
  plus per-energy diagnostics; see @ref calculate_formFactors_rpa.
*/
struct RPAFormFactors {
  //! Per core orbital (indexed as HF::HartreeFock::core()): the 13 factors,
  //! in the order and with the empty-if-not-requested convention of
  //! @ref calculate_formFactors_nk
  std::vector<std::array<LinAlg::Matrix<double>, 13>> K_nk{};
  //! Per energy deposit, worst over (rank, parity, q): RPA convergence
  //! achieved, iterations used, the |K - pi*D| internal-consistency
  //! deviation (ExternalField::TDHFcntm::KpiD_dev), and the K-matrix
  //! asymmetry. Zero for rows that were not dressed.
  std::vector<double> rpa_eps{}, rpa_its{}, KpiD_dev{}, Kbar_asym{};
};

/*!
  @brief Calculates the 13 form factors with core polarisation (RPA), for
  every core orbital at once.
  @details
  Same operator set, factor order, and (E, q) grids as
  @ref calculate_formFactors_nk, with every matrix element dressed by the
  continuum TDHF (RRPA) of ExternalField::TDHFcntm. The physical
  (unitarised) amplitude A of each open channel replaces the bare reduced
  matrix element,
  \f[ K(E,q) = (2k+1)\, x_a \sum_{\kappa_e} |A_{a\kappa_e}|^2 , \f]
  and the interference factors (X, Y, Z) use Re(A conj(A')) of the two
  operators in the same channel, whose relative phase is physical (see
  ExternalField::TDHFcntm::A_phys). The unitarised amplitudes need no HF
  continuum bra at all.

  The RPA couples all core orbitals, so one solve serves every hole
  orbital: the result is per orbital, not per call. Everything is computed
  at the outermost loop level it depends on:
  - core only: the RPA solver itself (one ExternalField::TDHFcntm instance,
    switched between operators, so its core-only tables are built once);
  - per E: the channel structure, and (standing-wave option only) the HF
    continuum bra of each ionised orbital;
  - per (E, rank, parity): the on-shell K matrix -- one field-free solve
    per open channel, shared by every operator of that rank and parity and
    by the whole q grid;
  - per (E, rank, operator, q): the driven solve. This is the dominant
    cost, and is parallelised over q (see RPAOptions::parallel_q).
  Energies run serially.

  @param vHF     Hartree-Fock potential (defines the core). HF states only.
  @param ec_min  Minimum continuum-electron energy to include (au).
  @param ec_max  Maximum continuum-electron energy to include (au).
  @param Egrid   Energy deposits E (au).
  @param qgrid   Momentum transfers q (au); one point if @p diagonal_Eq.
  @param diagonal_Eq  Fix q = E/c (absorption of a massless particle).
  @param low_q   Use the low-q forms of the operators.
  @param jK_tab  Spherical Bessel lookup table, ranks to at least Kmax.
  @param Kmin    Minimum multipole rank.
  @param Kmax    Maximum multipole rank.
  @param vectorQ, axialQ, scalarQ, pseudoscalarQ, spatialQ  Which factors,
                 as calculate_formFactors_nk.
  @param options  RPA options; see @ref RPAOptions.
  @param print   Print a per-energy progress and diagnostics line.
  @return See @ref RPAFormFactors.
  @note Rows above RPAOptions::E_max (and rows where no orbital is ionised)
  are left ZERO; the caller supplies the bare values there.
  @note The standing-wave option uses the orthogonalised V^{N-1} HF
  continuum bra of each orbital (hole_particle and force_orthog on), as
  TDHFcntm::dV_cntm requires; force_rescale does not apply.
*/
RPAFormFactors
calculate_formFactors_rpa(const HF::HartreeFock *vHF, double ec_min,
                          double ec_max, const std::vector<double> &Egrid,
                          const std::vector<double> &qgrid, bool diagonal_Eq,
                          bool low_q, const SphericalBessel::JL_table &jK_tab,
                          int Kmin, int Kmax, bool vectorQ, bool axialQ,
                          bool scalarQ, bool pseudoscalarQ, bool spatialQ,
                          const RPAOptions &options, bool print = true);

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