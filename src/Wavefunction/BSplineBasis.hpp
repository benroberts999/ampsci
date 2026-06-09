#pragma once
#include "IO/InputBlock.hpp"
#include "LinAlg/include.hpp"
#include "qip/String.hpp"
#include <memory>
#include <string>
#include <utility>
class DiracSpinor;
class Wavefunction;
class Grid;
namespace MBPT {
class CorrelationPotential;
}
namespace IO {
class InputBlock;
}

/*!
  @brief Constructs spinor/orbital basis using B-splines

  @details
  Uses Bsplines to form a set of B-spline orbitals (using method from [1]
  "Derevianko", [2] "Johnson", or [3] "Fischer"). 
  Diagonalises B-splines over Hamiltonian to produce a set of basis orbitals.

  * [1] K. Beloy, A. Derevianko, Comput. Phys. Commun. 179, 310 (2008).
  * [2] W. Johnson, S. Blundell, J. Sapirstein, Phys. Rev. A 37, 307 (1988).
  * [3] C. F. Fischer and F. A. Parpia, Phys. Lett. A 179, 198 (1993).
  * See also: 
    * Bachau et al., Reports Prog. Phys. 64, 1815 (2001).
    * V. M. Shabaev, I. I. Tupitsyn, V. A. Yerokhin, G. Plunien, and G. Soff, Phys. Rev. Lett. 93, 130405 (2004).

  By default, uses the Derevianko Dual Kintetic Balance (DKB) basis.

  If \f$\{|i\rangle\}\f$ are the set of \f$2N\f$ DKB spline orbitals (of a given
  angular symmetry), and
  \f[ 
    H_{ij} = \langle{S_i}|\hat H_{\rm HF}|{S_j}\rangle \,
    , \qquad S_{ij} = \langle{S_i|S_j}\rangle. 
  \f]
  The eigenvalue problem:
  \f[
    H_{ij}p_i = \epsilon S_{ij}p_i,
  \f]
  is solved, yielding \f$2N\f$ eigenvalues \f$\epsilon\f$ with corresponding
  eigenvectors \f$p\f$ (half of these are positive energy solutions, half are
  negative energy solutions).

  Note: form_basis() does not store the eigenvectors, instead, it expands the
  basis orbitals and stores them on the regular grid (coordinate space).
  i.e., for each eigenvalue, n, the corresponding basis orbital is:
  \f[
    |{n}\rangle = \sum_i^{2N} p_i |i\rangle\,
  \f]
*/
namespace SplineBasis {

/*!
  @brief B-spline construction method: Derevianko/Reno (DKB), Johnson/NotreDame, or Fischer.
  @details
  W. R. Johnson, S. A. Blundell, J. Sapirstein, Phys. Rev. A 37, 307 (1988).
*/
enum class SplineType { Derevianko, Johnson, Fischer };

//! Returns the canonical name string for a SplineType.
inline std::string_view parseSplineType(SplineType type) {
  switch (type) {
  case SplineType::Johnson:
    return "Johnson";
  case SplineType::Fischer:
    return "Fischer";
  case SplineType::Derevianko:
    return "Derevianko";
  }
  return "UnknownSplineType";
}

/*!
  @brief Parses a string to SplineType (case-insensitive).
  @details
  - Derevianko (default, aliases: reno): Dual Kinetic Balance method; K. Beloy, A. Derevianko, Comput. Phys. Commun. 179, 310 (2008).
  - Johnson (aliases: nd, notredame, notre-dame): Notre-Dame method with explicit boundary conditions; W. R. Johnson, S. A. Blundell, J. Sapirstein, Phys. Rev. A 37, 307 (1988).
  - Fischer (aliases: vanderbilt): rmax boundary only (no r=0 conditions); C. F. Fischer and F. A. Parpia, Phys. Lett. A 179, 198 (1993).
*/
inline auto parseSplineType(std::string_view type) {
  using qip::ci_compare;
  if (ci_compare(type, "Johnson") || ci_compare(type, "nd") ||
      ci_compare(type, "NotreDame") || ci_compare(type, "Notre-Dame"))
    return SplineType::Johnson;
  if (ci_compare(type, "Fischer") || ci_compare(type, "Vanderbilt"))
    return SplineType::Fischer;
  if (ci_compare(type, "Derevianko") || ci_compare(type, "Reno") ||
      type.empty())
    return SplineType::Derevianko;
  std::cout << "Warning: unknown SplineType '" << type
            << "'. Defaulting to Derevianko.\n";
  return SplineType::Derevianko;
}

//! Input parameters for B-spline basis construction.
struct Parameters {
  Parameters() {}
  Parameters(const std::string &states, std::size_t n, std::size_t k, double r0,
             double reps, double rmax, const std::string &positron = "",
             SplineType itype = SplineType::Derevianko,
             bool in_orthogonalise = false, bool in_verbose = true);
  Parameters(IO::InputBlock input);

  std::string states{};
  std::size_t n{}, k{};
  double r0{}, reps{}, rmax{};
  std::string positron{};
  SplineType type{SplineType::Derevianko};
  bool orthogonalise{false};
  bool verbose{true};
};

//!
/*!
  @brief Forms and returns the basis orbitals (expanded in terms of splines).
  @details
  - states_str = which states to keep e.g., "25spd10f" (up to n=25 for
  s,p,d-states, and up to n=10 for f states)
  - n_spl: Number of splines (nb: underlying spline set is larger, see [1])
  - k_spl: k order of the B-splines. NB: must have \f$k\geq l_{\rm max}+3\f$ [1]
  - r0_spl: first internal knot
  - r0_eps: sets r0_spl as r where relative core density larger than
 r0_eps (updates r0 for each l). Typically ~1.0e-8. Set to zero to use r0_spl.
  - rmax_spl: last internal knot (basis orbitals only non-zero before this
  point)
  - wf: Wavefunction object: needed to form Hartree-Fock Hamiltonian
  - positronQ: =true will keep negative energy states (have -ve principal
  quantum number, are appended to end of the basis std::vector). If false,
  discards them.

Note: This function calls the below functions, they rarely need to be called
explicitely, unless you are trying to do something different to usual.
*/
std::vector<DiracSpinor> form_basis(const Parameters &params,
                                    const Wavefunction &wf,
                                    const bool correlationsQ = false);

//! Compares basis to reference orbitals; checks normality, orthogonality, and energies.
double check(const std::vector<DiracSpinor> &basis,
             const std::vector<DiracSpinor> &orbs, bool print_warning = true);

//! Type-dependent spline grid parameters for a given kappa.
struct SplineParams {
  //! first included spline index
  std::size_t imin;
  //! total number of splines (n_spl >= N, number of basis functions)
  std::size_t n_spl;
  //! one-past-last included spline index (exclusive upper bound)
  std::size_t imax;
  //! kinetic balance prefactor (1 for DKB, 0 for Johnson)
  double lambda_DKB;
};
//! Returns SplineParams for the given type, kappa, and number of states.
SplineParams spline_params(SplineType type, int kappa, std::size_t n_states);

//! Forms the underlying spline basis (which is not kept)
std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>> form_spline_basis(
  const int kappa, const std::size_t n_states, const std::size_t k_spl,
  const double r0_spl, const double rmax_spl, std::shared_ptr<const Grid> rgrid,
  const double alpha, SplineType itype = SplineType::Derevianko);

//! Calculates and returns the Hamiltonian matrix H_ij and overlap matrix S_ij.
std::pair<LinAlg::Matrix<double>, LinAlg::Matrix<double>>
fill_Hamiltonian_matrix(const std::vector<DiracSpinor> &spl_basis,
                        const std::vector<DiracSpinor> &d_basis,
                        const Wavefunction &wf,
                        const bool correlationsQ = false,
                        SplineType itype = SplineType::Derevianko);

/*!
  @brief Applies Johnson/Notre-Dame boundary conditions to the Hamiltonian matrix.
  @details
  Adds r=0 and r=rmax boundary terms to enforce correct asymptotic behaviour.
  W. R. Johnson, S. A. Blundell, J. Sapirstein, Phys. Rev. A 37, 307 (1988).
*/
void add_JohnsonBoundary(LinAlg::Matrix<double> *Aij, const int kappa,
                         const double alpha);

/*!
  @brief Applies Fischer/Vanderbilt rmax boundary: forces f(rmax)=g(rmax).
  @details
  Subset of add_JohnsonBoundary (rmax terms only; no r=0 conditions).
*/
void add_FischerBoundary(LinAlg::Matrix<double> *Aij, const double alpha);

//! Expands basis orbitals in terms of spline orbitals by diagonalising Hamiltonian.
void expand_basis_orbitals(std::vector<DiracSpinor> *basis,
                           std::vector<DiracSpinor> *basis_positron,
                           const std::vector<DiracSpinor> &spl_basis,
                           const int kappa, const int max_n, int max_n_positron,
                           const LinAlg::Vector<double> &e_values,
                           const LinAlg::Matrix<double> &e_vectors,
                           const Wavefunction &wf);

//! TKR sum rule (basis test); should =0 (must include -ve energy states)
std::vector<double> sumrule_TKR(const std::vector<DiracSpinor> &basis,
                                const std::vector<double> &r,
                                bool print = false);

//! Drake-Gordon sum rule (basis test); should =0 (must incl -ve energy states)
std::vector<double> sumrule_DG(int nDG, const std::vector<DiracSpinor> &basis,
                               const Grid &gr, double alpha, bool print);

//! Measures radial completeness of the basis for orbital Fv.
std::pair<double, double> r_completeness(const DiracSpinor &Fv,
                                         const std::vector<DiracSpinor> &basis,
                                         const Grid &gr, bool print = false);

//! Prints the title and column headers for r_completeness output.
inline void r_completeness_header() {
  std::printf("Completeness (radial sum rules):\n");
  std::printf("  <a|a>    = sum_n <a|r|n><n|1/r|a>  [= 1]\n");
  std::printf("  <a|r2|a> = sum_n <a|r|n><n|r|a>    [= <a|r^2|a>_HF]\n");
  std::printf("%-4s  %7s [%7s] %6s | %10s [%10s] %6s\n", "St", "sum1", "1",
              "eps", "sum_r2", "<r2>_HF", "eps");
}
} // namespace SplineBasis
