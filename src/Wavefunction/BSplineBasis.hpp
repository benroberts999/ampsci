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
    * V. M. Shabaev, I. I. Tupitsyn, V. A. Yerokhin, G. Plunien, and G. Soff,
      Phys. Rev. Lett. 93, 130405 (2004).

  By default, uses the Derevianko Dual Kinetic Balance (DKB) basis.

  If \f$\{|i\rangle\}\f$ are the set of \f$2N\f$ DKB spline orbitals (of a given
  angular symmetry), and
  \f[
    H_{ij} = \langle{S_i}|\hat H_{\rm HF}|{S_j}\rangle \,
    , \qquad S_{ij} = \langle{S_i|S_j}\rangle.
  \f]
  The eigenvalue problem:
  \f[
    H_{ij}p_j = \epsilon S_{ij}p_j,
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
  @brief B-spline construction method: Derevianko (DKB), Johnson, or Fischer.

  @details
  - Derevianko (Reno): K. Beloy, A. Derevianko, Comput. Phys. Commun. 179,
    310 (2008).
  - Johnson (Notre Dame): W. R. Johnson, S. A. Blundell, J. Sapirstein,
    Phys. Rev. A 37, 307 (1988).
  - Fischer (Vanderbilt): C. F. Fischer and F. A. Parpia, Phys. Lett. A 179,
    198 (1993).

  See parseSplineType for a description of each.
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
  - Derevianko (default, aliases: reno): Dual Kinetic Balance method;
    K. Beloy, A. Derevianko, Comput. Phys. Commun. 179, 310 (2008).
  - Johnson (aliases: nd, notredame, notre-dame): Notre-Dame method with
    explicit boundary conditions; W. R. Johnson, S. A. Blundell,
    J. Sapirstein, Phys. Rev. A 37, 307 (1988).
  - Fischer (aliases: vanderbilt): rmax boundary only (no r=0 conditions);
    C. F. Fischer and F. A. Parpia, Phys. Lett. A 179, 198 (1993).

  @return The SplineType; Derevianko if @p type is empty.

  @warning An unrecognised @p type prints a warning and returns Derevianko.
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

//! Input parameters for B-spline basis construction (see form_basis).
struct Parameters {
  //! Default: all parameters empty/zero
  Parameters() {}
  //! Sets each parameter directly (see form_basis)
  Parameters(const std::string &states, std::size_t n, std::size_t k, double r0,
             double reps, double rmax, const std::string &positron = "",
             SplineType itype = SplineType::Derevianko,
             bool in_orthogonalise = false, bool in_verbose = true);
  /*!
    @brief Reads parameters from an input block.

    @details
    Input keys [defaults]: states [""], number (n) [30], order (k) [7],
    r0 [1.0e-4], r0_eps (reps) [0.0], rmax [40.0], positron [""],
    type [Derevianko], orthogonalise [false], verbose [true].
  */
  Parameters(IO::InputBlock input);

  //! Electron states to keep, e.g., "25spd10f"
  std::string states{};
  //! Number of splines (n), and order of the B-splines (k)
  std::size_t n{}, k{};
  //! First internal knot (r0), core-density cut-off (reps), last knot (rmax)
  double r0{}, reps{}, rmax{};
  //! Negative-energy states to keep, e.g., "9spdf" (empty: none)
  std::string positron{};
  //! B-spline construction method
  SplineType type{SplineType::Derevianko};
  //! Force orthogonal to core (spectrum: to valence); done by Wavefunction
  bool orthogonalise{false};
  //! Print construction details to screen
  bool verbose{true};
};

/*!
  @brief Forms and returns the basis orbitals (expanded in terms of splines).

  @details
  Members of @p params:
  - states = which states to keep e.g., "25spd10f" (up to n=25 for
    s,p,d-states, and up to n=10 for f states)
  - n: Number of splines (nb: underlying spline set is larger, see [1])
  - k: k order of the B-splines. NB: must have \f$k\geq l_{\rm max}+3\f$ [1]
  - r0: first internal knot
  - reps: sets r0 as r where relative core density larger than reps, if this
    is larger than r0 (updates r0 for each l; for l above the core, uses
    the largest core l). Typically ~1.0e-8. Set to zero to use r0.
  - rmax: last knot (basis orbitals only non-zero before this point)
  - positron: basis string giving the number of negative-energy states to
    keep for each kappa (e.g., "9spdf" keeps 9 of each). Those closest to
    -2mc^2 are kept. They are labelled with n continuing on from the
    positive-energy solutions (so {n,kappa} stays unique), are flagged
    (DiracSpinor::negativeEnergyStateQ), and are appended to the end of the
    basis std::vector. Empty string discards them.
  - type: B-spline construction method (see SplineType)
  - verbose: print construction details to screen
  - orthogonalise: not used here (see Wavefunction::formBasis)

  @param params         Basis parameters (see above)
  @param wf             Wavefunction object: needed to form Hartree-Fock
                        Hamiltonian
  @param correlationsQ  If true, includes the correlation potential (Sigma)
                        of @p wf in the Hamiltonian (if it exists)

  @return Basis orbitals: electron states for each kappa (in the order of
          @p params.states, increasing n), followed by any negative-energy
          states.

  @note This function calls the below functions, they rarely need to be
        called explicitly, unless you are trying to do something different
        to usual.
*/
std::vector<DiracSpinor> form_basis(const Parameters &params,
                                    const Wavefunction &wf,
                                    const bool correlationsQ = false);

/*!
  @brief Checks basis against orbitals: normality, orthogonality, and energies.

  @details
  For each orbital in @p orbs with a matching basis state (same n and kappa),
  finds the worst \f$ ||\langle a|b\rangle| - 1| \f$ and relative energy
  difference, \f$ |(\varepsilon_b - \varepsilon_a)/\varepsilon_b| \f$. Also
  finds the worst orthonormality between @p orbs and @p basis states of the
  same kappa (see DiracSpinor::check_ortho). Prints each, and lists basis
  states with opposite sign to the matching orbital.

  @param basis          Basis orbitals
  @param orbs           Reference orbitals (e.g., core or valence)
  @param print_warning  Flag printed results worse than 1e-3

  @return Largest of the three deviations; 0 (prints nothing) if @p orbs is
          empty.
*/
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
  //! kinetic balance prefactor (1 for DKB and Fischer, 0 for Johnson)
  double lambda_DKB;
};

/*!
  @brief Returns SplineParams for the given type, kappa, and number of states.

  @details
  Of the n_spl underlying B-splines, @p n_states are kept (imin to imax-1):
  - Johnson    : all; no kinetic balance
  - Fischer    : drops the first; kinetic balance
  - Derevianko : drops the first |kappa| and the last; kinetic balance
*/
SplineParams spline_params(SplineType type, int kappa, std::size_t n_states);

/*!
  @brief Forms the underlying spline basis (which is not kept).

  @details
  Forms 2N spline spinors for given kappa (N = @p n_states), from the
  B-splines \f$B_i\f$ selected by spline_params: N "f-like" and N "g-like",
  \f[
    S_i = \begin{pmatrix} B_i \\
      \frac{\lambda\alpha}{2}\left(B_i' + \frac{\kappa}{r}B_i\right)
    \end{pmatrix},
    \qquad
    S_{N+i} = \begin{pmatrix}
      \frac{\lambda\alpha}{2}\left(B_i' - \frac{\kappa}{r}B_i\right) \\ B_i
    \end{pmatrix},
  \f]
  where \f$\lambda\f$ is the kinetic balance prefactor (SplineParams). Knots
  are distributed log-linearly between @p r0_spl and @p rmax_spl.

  @param kappa     Dirac quantum number
  @param n_states  Number of splines, N, of each type (f-like and g-like)
  @param k_spl     Order of the B-splines
  @param r0_spl    First internal knot
  @param rmax_spl  Last knot (splines are zero beyond this point)
  @param rgrid     Radial grid the splines are stored on
  @param alpha     Fine-structure constant
  @param itype     B-spline construction method

  @return Pair {splines, derivatives}: the 2N spline spinors, and their
          radial derivatives (d/dr of both f and g).
*/
std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>> form_spline_basis(
  const int kappa, const std::size_t n_states, const std::size_t k_spl,
  const double r0_spl, const double rmax_spl, std::shared_ptr<const Grid> rgrid,
  const double alpha, SplineType itype = SplineType::Derevianko);

/*!
  @brief Calculates and returns the Hamiltonian matrix H_ij and overlap S_ij.

  @details
  \f$ H_{ij} = \langle S_i|\hat H|S_j\rangle \f$ and
  \f$ S_{ij} = \langle S_i|S_j\rangle \f$, over the spline basis. The
  Hamiltonian is that of @p wf: the local part (see Wavefunction::H0ab),
  plus Hartree-Fock exchange (if non-local), Breit (if included in @p wf),
  and Sigma (if @p correlationsQ and @p wf has Sigma). For Johnson and
  Fischer types, also adds the boundary terms (see add_JohnsonBoundary and
  add_FischerBoundary).

  @param spl_basis      Spline basis (see form_spline_basis)
  @param d_basis        Radial derivatives of the spline basis
  @param wf             Wavefunction: provides the Hamiltonian
  @param correlationsQ  Include the correlation potential (Sigma) of @p wf
  @param itype          B-spline construction method (sets boundary terms)

  @return Pair {H_ij, S_ij}.
*/
std::pair<LinAlg::Matrix<double>, LinAlg::Matrix<double>>
fill_Hamiltonian_matrix(const std::vector<DiracSpinor> &spl_basis,
                        const std::vector<DiracSpinor> &d_basis,
                        const Wavefunction &wf,
                        const bool correlationsQ = false,
                        SplineType itype = SplineType::Derevianko);

/*!
  @brief Applies Johnson/Notre-Dame boundary conditions to Hamiltonian matrix.

  @details
  Adds r=0 and r=rmax boundary terms to enforce correct asymptotic behaviour.
  W. R. Johnson, S. A. Blundell, J. Sapirstein, Phys. Rev. A 37, 307 (1988).

  @param Aij    Hamiltonian matrix, 2N x 2N (N f-like splines, then N
                g-like); updated in place
  @param kappa  Dirac quantum number
  @param alpha  Fine-structure constant
*/
void add_JohnsonBoundary(LinAlg::Matrix<double> *Aij, const int kappa,
                         const double alpha);

/*!
  @brief Applies Fischer/Vanderbilt rmax boundary: forces f(rmax)=g(rmax).

  @details
  Subset of add_JohnsonBoundary (rmax terms only; no r=0 conditions).

  @param Aij    Hamiltonian matrix, 2N x 2N (N f-like splines, then N
                g-like); updated in place
  @param alpha  Fine-structure constant
*/
void add_FischerBoundary(LinAlg::Matrix<double> *Aij, const double alpha);

/*!
  @brief Expands basis orbitals in terms of spline orbitals.

  @details
  Basis orbitals are formed from the eigenvectors of the Hamiltonian.
  Keeps the electron states up to n = @p max_n (labelled n = l+1, l+2, ...,
  from the bottom of the positive-energy branch up), and the
  @p max_n_positron negative-energy states closest to the gap (labelled with
  n continuing on from the positive-energy solutions, and flagged).

  Each orbital is normalised, with sign chosen so that f is positive at
  small r.

  @param basis           Electron states are appended to this
  @param basis_positron  Negative-energy states are appended to this
  @param spl_basis       Spline basis the eigenvectors are expressed in
  @param kappa           Dirac quantum number
  @param max_n           Maximum n of electron states kept
  @param max_n_positron  Number of negative-energy states kept
  @param e_values        Eigenvalues (ascending)
  @param e_vectors       Eigenvectors: e_vectors[i] corresponds to e_values[i]
  @param wf              Wavefunction: provides grid and alpha

  @note States whose norm differs from 1 by more than 1e-4 are treated as
        spurious: a warning is printed, and the orbital is set to zero (but
        is still kept).
*/
void expand_basis_orbitals(std::vector<DiracSpinor> *basis,
                           std::vector<DiracSpinor> *basis_positron,
                           const std::vector<DiracSpinor> &spl_basis,
                           const int kappa, const int max_n, int max_n_positron,
                           const LinAlg::Vector<double> &e_values,
                           const LinAlg::Matrix<double> &e_vectors,
                           const Wavefunction &wf);

/*!
  @brief TKR sum rule (basis test); should be zero (needs -ve energy states).

  @details
  Thomas-Kuhn-Reiche sum rule. With a the first state in @p basis, for each
  l up to the maximum l in @p basis:
  \f[
    \sum_n \frac{w_n}{2l+1} (\varepsilon_n - \varepsilon_a)
    \langle a|r|n\rangle^2 = 0,
  \f]
  summed over basis states with \f$\kappa_n = l\f$ (\f$w_n = l\f$) and
  \f$\kappa_n = -l-1\f$ (\f$w_n = l+1\f$), where
  \f$\langle a|r|n\rangle\f$ is the radial integral.

  @param basis  Basis orbitals (including negative-energy states)
  @param r      Radial grid points
  @param print  Print electron and negative-energy contributions, for each l

  @return The sum, for each l; empty if @p basis is empty.
*/
std::vector<double> sumrule_TKR(const std::vector<DiracSpinor> &basis,
                                const std::vector<double> &r,
                                bool print = false);

/*!
  @brief Drake-Goldman sum rules (basis test); should be zero.

  @details
  Must include negative-energy states in @p basis.
  For each kappa in @p basis (a: the first basis state with that kappa):
  \f[
    \frac{1}{2j_a+1}\sum_n (\varepsilon_n - \varepsilon_a)^{m}
    |\langle a||r||n\rangle|^2
  \f]
  where m = @p nDG, compared with its exact value:
  \f$\langle a|r^2|a\rangle\f$ (m=0) or 0 (m=1). For m=2, the sum is scaled
  by \f$\alpha^2/3\f$, and compared with 1. Kappas with l equal to the
  maximum l in @p basis are skipped (no l+1 states to sum over).

  @param nDG    Power of the energy difference, m: 0, 1, or 2
  @param basis  Basis orbitals (including negative-energy states)
  @param gr     Radial grid
  @param alpha  Fine-structure constant
  @param print  Print sum, exact value, and difference, for each kappa

  @return Sum minus exact value, for each kappa included.

  @warning Every kappa up to the largest in @p basis must be present (e.g.,
           not "20s5d", which has no p states): a missing kappa dereferences
           an end iterator.
*/
std::vector<double> sumrule_DG(int nDG, const std::vector<DiracSpinor> &basis,
                               const Grid &gr, double alpha, bool print);

/*!
  @brief Measures radial completeness of the basis for orbital Fv.

  @details
  Sums over basis states n with the same kappa as @p Fv:
  \f[
    \sum_n \langle v|r|n\rangle\langle n|r^{-1}|v\rangle = 1, \qquad
    \sum_n \langle v|r|n\rangle\langle n|r|v\rangle = \langle v|r^2|v\rangle,
  \f]
  where the matrix elements are radial integrals.

  @param Fv     Orbital to test
  @param basis  Basis orbitals
  @param gr     Radial grid
  @param print  Print one line of results (see r_completeness_header)

  @return Relative deviation of each sum from its exact value.
*/
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
