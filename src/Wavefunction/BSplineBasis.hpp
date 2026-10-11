#pragma once
#include "IO/InputBlock.hpp"
#include "LinAlg/include.hpp"
#include "qip/String.hpp"
#include <memory>
#include <string>
#include <utility>
#include <vector>
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
  "Derevianko", [2] "Shabaev", [3] "Johnson", or [4] "Fischer", or the
  "Hybrid" of [1] and [3]; see SplineType).
  Diagonalises B-splines over Hamiltonian to produce a set of basis orbitals.

  * [1] K. Beloy, A. Derevianko, Comput. Phys. Commun. 179, 310 (2008).
  * [2] V. M. Shabaev, I. I. Tupitsyn, V. A. Yerokhin, G. Plunien, and G. Soff,
    Phys. Rev. Lett. 93, 130405 (2004).
  * [3] W. Johnson, S. Blundell, J. Sapirstein, Phys. Rev. A 37, 307 (1988).
  * [4] C. F. Fischer and F. A. Parpia, Phys. Lett. A 179, 198 (1993).
  * See also:
    * Bachau et al., Reports Prog. Phys. 64, 1815 (2001).

  By default, uses the Derevianko Dual Kinetic Balance (DKB) basis.

  If \f$\{|i\rangle\}\f$ are the set of spline orbitals (of a given angular
  symmetry; \f$2N\f$ of them, except Fischer with \f$\kappa>0\f$, see
  spline_params), and
  \f[
    H_{ij} = \langle{S_i}|\hat H_{\rm HF}|{S_j}\rangle \,
    , \qquad S_{ij} = \langle{S_i|S_j}\rangle.
  \f]
  The eigenvalue problem:
  \f[
    H_{ij}p_j = \epsilon S_{ij}p_j,
  \f]
  is solved, yielding one eigenvalue \f$\epsilon\f$ (with eigenvector
  \f$p\f$) for each spline orbital: roughly half are positive-energy
  solutions, and half negative-energy solutions.

  Note: form_basis() does not store the eigenvectors, instead, it expands the
  basis orbitals and stores them on the regular grid (coordinate space).
  i.e., for each eigenvalue, n, the corresponding basis orbital is:
  \f[
    |{n}\rangle = \sum_i p_i |i\rangle\,
  \f]
*/
namespace SplineBasis {

/*!
  @brief B-spline construction method (see details).

  @details
  - Derevianko (Reno): dual kinetic balance (DKB) splines; the first |kappa|
    and the last spline are dropped, so the B-spline parts vanish at both
    ends of the cavity (the kinetic-balance parts need not).
    K. Beloy, A. Derevianko, Comput. Phys. Commun. 179, 310 (2008).
  - Shabaev: DKB splines; the first and last spline are dropped, for all
    kappa, with no boundary terms (as published: dropping more at r = 0 for
    |kappa| >= 2 was found to make practically no difference).
    V. M. Shabaev, I. I. Tupitsyn, V. A. Yerokhin, G. Plunien, and G. Soff,
    Phys. Rev. Lett. 93, 130405 (2004).
  - Johnson (Notre Dame): plain splines for f and g (no kinetic balance),
    all kept; the boundary conditions f(0) = 0 and f(R) = -g(R) (MIT bag)
    are imposed through boundary terms in the Hamiltonian matrix.
    W. R. Johnson, S. A. Blundell, J. Sapirstein, Phys. Rev. A 37, 307
    (1988).
  - Fischer: plain splines (no kinetic balance); the first and
    last spline are dropped (f = g = 0 at both ends), with no boundary
    terms. For kappa > 0, the second spline is also dropped from f
    (f'(0) = 0), which removes the spurious state found by Fischer and
    Parpia for kappa = 1.
    C. F. Fischer and F. A. Parpia, Phys. Lett. A 179, 198 (1993).
  - Hybrid: DKB splines; the first |kappa| are dropped (as for
    Derevianko), but the last is kept, and the Johnson rmax (MIT bag)
    boundary term f(R) = -g(R) is imposed (see add_rmaxBoundary). Differs
    from Derevianko only in the boundary condition at rmax.

  See parseSplineType for the input strings.

  @warning Fischer: bound states are accurate, but the pseudo-continuum
  (cavity) states are not reliable: they disagree with those of the other
  types, and for kappa > 0 some lie below the corresponding free-atom
  energies. May not be suitable for MBPT; intended for testing only.
*/
enum class SplineType { Derevianko, Shabaev, Johnson, Fischer, Hybrid };

//! Returns the canonical name string for a SplineType.
inline std::string_view parseSplineType(SplineType type) {
  switch (type) {
  case SplineType::Johnson:
    return "Johnson";
  case SplineType::Fischer:
    return "Fischer";
  case SplineType::Hybrid:
    return "Hybrid";
  case SplineType::Shabaev:
    return "Shabaev";
  case SplineType::Derevianko:
    return "Derevianko";
  }
  return "UnknownSplineType";
}

/*!
  @brief Parses a string to SplineType (case-insensitive).

  @details
  - Derevianko (default; aliases: reno, dkb)
  - Shabaev
  - Johnson (aliases: nd, notredame, notre-dame)
  - Fischer
  - Hybrid

  See SplineType for a description of each.

  @return The SplineType; Derevianko if @p type is empty.

  @warning An unrecognised @p type prints a warning and returns Derevianko.
*/
inline auto parseSplineType(std::string_view type) {
  using qip::ci_compare;
  if (ci_compare(type, "Johnson") || ci_compare(type, "nd") ||
      ci_compare(type, "NotreDame") || ci_compare(type, "Notre-Dame"))
    return SplineType::Johnson;
  if (ci_compare(type, "Fischer"))
    return SplineType::Fischer;
  if (ci_compare(type, "Hybrid"))
    return SplineType::Hybrid;
  if (ci_compare(type, "Shabaev"))
    return SplineType::Shabaev;
  if (ci_compare(type, "Derevianko") || ci_compare(type, "Reno") ||
      ci_compare(type, "DKB") || type.empty())
    return SplineType::Derevianko;
  std::cout << "Warning: unknown SplineType '" << type
            << "'. Defaulting to Derevianko.\n";
  return SplineType::Derevianko;
}

/*!
  @brief Default first internal knot, r0, for each l, for the given type.

  @details
  Element l is the first knot for orbital angular momentum l; the last
  element is used for all higher l (see form_basis).
  - Dual-kinetic-balance types (Derevianko, Shabaev, Hybrid): 1e-4 for all l
  - Types without kinetic balance (Johnson, Fischer): 3e-4 for s, 1e-4 for
    p, 1e-3 for d and f, 1e-2 for higher l. As in Johnson et al. (Table
    III), the first knot is moved out for higher l, since these orbitals are
    negligible close to the nucleus, which improves the knot distribution.
    The s and p values give the most accurate hyperfine constants and
    Drake-Goldman sum rules for Cs (see the basis tests); in particular, a
    first s knot at 1e-4 degrades the s-state hyperfine constants.
*/
std::vector<double> default_r0(SplineType type);

//! Input parameters for B-spline basis construction (see form_basis).
struct Parameters {
  //! Default: all parameters empty/zero
  Parameters() {}
  //! Sets each parameter directly (see form_basis); empty r0: default_r0
  Parameters(const std::string &states, std::size_t n, std::size_t k,
             std::vector<double> r0, double rmax,
             const std::string &positron = "",
             SplineType itype = SplineType::Derevianko,
             bool in_orthogonalise = false, bool in_verbose = true);
  /*!
    @brief Reads parameters from an input block.

    @details
    Input keys [defaults]: states [""], number (n) [30], order (k) [7],
    r0 [default_r0(type); comma-separated list, one for each l], rmax [40.0],
    positron [""], type [Derevianko], orthogonalise [false], verbose [true].
  */
  Parameters(IO::InputBlock input);

  //! Electron states to keep, e.g., "25spd10f"
  std::string states{};
  //! Number of splines (n), and order of the B-splines (k)
  std::size_t n{}, k{};
  //! First internal knot for each l (last entry for higher l); empty: default
  std::vector<double> r0{};
  //! Last knot (rmax)
  double rmax{};
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
  - n: Number of splines (the underlying spline set may be larger; see
    spline_params)
  - k: k order of the B-splines. NB: must have \f$k\geq l_{\rm max}+3\f$ [1]
  - r0: first internal knot, for each l: r0[l] is used for orbital angular
    momentum l, and the last entry for all higher l (so a single value
    applies to all l). Empty: default_r0(type). Values below the first grid
    point are replaced by it.
  - rmax: last knot (basis orbitals only non-zero before this point)
  - positron: basis string giving the number of negative-energy states to
    keep for each kappa (e.g., "9spdf" keeps 9 of each). Those closest to
    -2mc^2 are kept. They are labelled with n continuing on from the
    positive-energy solutions (so {n,kappa} stays unique), are flagged
    (DiracSpinor::negativeEnergyStateQ), and are appended to the end of the
    basis std::vector. Empty string discards them. Must be a subset of states.
    A warning is also printed if any label is too large for DiracSpinor::nk_index (n > 256;
    e.g., many splines with many negative-energy states kept).
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

//! Type-dependent spline selection for a given kappa (see spline_params).
struct SplineParams {
  //! total number of underlying B-splines
  std::size_t n_spl;
  //! f-like splines kept: indices [imin_f, imax_f)
  std::size_t imin_f, imax_f;
  //! g-like splines kept: indices [imin_g, imax_g)
  std::size_t imin_g, imax_g;
  //! kinetic balance prefactor (1 for DKB types, 0 for Johnson, Fischer)
  double lambda_DKB;
  //! number of f-like splines kept
  std::size_t num_f() const { return imax_f - imin_f; }
  //! number of g-like splines kept
  std::size_t num_g() const { return imax_g - imin_g; }
};

/*!
  @brief Returns SplineParams for the given type, kappa, and number of states.

  @details
  Of the n_spl underlying B-splines, @p n_states are kept for each of f and
  g (indices [imin_f, imax_f) and [imin_g, imax_g)):
  - Johnson    : all; no kinetic balance
  - Fischer    : drops the first and the last; no kinetic balance. For
                 kappa > 0, the second is also dropped from f (f'(0) = 0),
                 which removes the spurious kappa = 1 state found by Fischer
                 and Parpia; f then has @p n_states - 1 splines
  - Hybrid     : drops the first |kappa|; kinetic balance
  - Shabaev    : drops the first and the last; kinetic balance
  - Derevianko : drops the first |kappa| and the last; kinetic balance
*/
SplineParams spline_params(SplineType type, int kappa, std::size_t n_states);

//! Values at r of the two spline spinors built from one B-spline B_i
struct SplineSpinorValues {
  //! f-like spinor (f = B_i): f, g, and their radial derivatives
  double f_f, g_f, df_f, dg_f;
  //! g-like spinor (g = B_i): f, g, and their radial derivatives
  double f_g, g_g, df_g, dg_g;
};

/*!
  @brief Values at r of the f-like and g-like spline spinors built from B_i.

  @details
  See form_spline_basis: the f-like spinor has f = B_i and
  g = kb (B_i' + kappa B_i/r); the g-like spinor has g = B_i and
  f = kb (B_i' - kappa B_i/r), where kb = lambda alpha/2 is the kinetic
  balance factor (zero without kinetic balance).

  @param B      B_i(r)
  @param dB     B_i'(r)
  @param d2B    B_i''(r)
  @param kappa  Dirac quantum number
  @param r      Radius (> 0)
  @param kb     Kinetic balance factor, lambda alpha/2
*/
SplineSpinorValues spline_spinor_values(double B, double dB, double d2B,
                                        int kappa, double r, double kb);

//! The spline spinors for one kappa (see form_spline_basis)
struct SplineSpinors {
  //! Spline spinors on the radial grid: f-like, then g-like
  std::vector<DiracSpinor> basis{};
  //! Radial derivatives (d/dr of both f and g) of the spline spinors
  std::vector<DiracSpinor> d_basis{};
  //! f and g of each spline spinor exactly at rmax (for boundary terms)
  std::vector<double> f_rmax{}, g_rmax{};
};

/*!
  @brief Forms the underlying spline basis (which is not kept).

  @details
  Forms the spline spinors for given kappa from the B-splines \f$B_i\f$
  selected by spline_params: \f$N_f\f$ "f-like" followed by \f$N_g\f$
  "g-like" (both equal to @p n_states, except Fischer with kappa > 0),
  \f[
    S_i = \begin{pmatrix} B_i \\
      \frac{\lambda\alpha}{2}\left(B_i' + \frac{\kappa}{r}B_i\right)
    \end{pmatrix},
    \qquad
    S_{N_f+i} = \begin{pmatrix}
      \frac{\lambda\alpha}{2}\left(B_i' - \frac{\kappa}{r}B_i\right) \\ B_i
    \end{pmatrix},
  \f]
  where \f$\lambda\f$ is the kinetic balance prefactor (SplineParams). Knots
  are distributed log-linearly between @p r0_spl and @p rmax_spl, with the
  end points (0 and rmax) repeated k times: only the first spline is
  non-zero at r = 0, and only the last at r = rmax.

  @param kappa     Dirac quantum number
  @param n_states  Number of splines, N, of each type (f-like and g-like;
                   N-1 f-like for Fischer with kappa > 0)
  @param k_spl     Order of the B-splines
  @param r0_spl    First internal knot
  @param rmax_spl  Last knot (splines are zero beyond this point)
  @param rgrid     Radial grid the splines are stored on
  @param alpha     Fine-structure constant
  @param itype     B-spline construction method

  @return The N_f + N_g spline spinors, their radial derivatives (d/dr of
          both f and g), and the values of f and g of each exactly at rmax
          (which is generally between grid points).
*/
SplineSpinors form_spline_basis(const int kappa, const std::size_t n_states,
                                const std::size_t k_spl, const double r0_spl,
                                const double rmax_spl,
                                std::shared_ptr<const Grid> rgrid,
                                const double alpha,
                                SplineType itype = SplineType::Derevianko);

/*!
  @brief Calculates and returns the Hamiltonian matrix H_ij and overlap S_ij.

  @details
  \f$ H_{ij} = \langle S_i|\hat H|S_j\rangle \f$ and
  \f$ S_{ij} = \langle S_i|S_j\rangle \f$, over the spline basis. The
  Hamiltonian is that of @p wf: the local part (see Wavefunction::H0ab),
  plus Hartree-Fock exchange (if non-local), Breit (if included in @p wf),
  and Sigma (if @p correlationsQ and @p wf has Sigma). For Johnson and
  Hybrid types, also adds the boundary terms (see add_JohnsonBoundary and
  add_rmaxBoundary).

  @param splines        Spline spinors (see form_spline_basis)
  @param wf             Wavefunction: provides the Hamiltonian
  @param correlationsQ  Include the correlation potential (Sigma) of @p wf
  @param itype          B-spline construction method (sets boundary terms)

  @return Pair {H_ij, S_ij}.
*/
std::pair<LinAlg::Matrix<double>, LinAlg::Matrix<double>>
fill_Hamiltonian_matrix(const SplineSpinors &splines, const Wavefunction &wf,
                        const bool correlationsQ = false,
                        SplineType itype = SplineType::Derevianko);

/*!
  @brief Applies Johnson/Notre-Dame boundary conditions to Hamiltonian matrix.

  @details
  Adds the boundary term of the action, S' [Eq. (10) of Johnson et al.],
  together with the surface term \f$\frac{c}{2}[fg]_0^R\f$ by which the
  symmetric kinetic matrix (Wavefunction::H0ab) differs from the
  antisymmetrised form of Johnson et al. [Eq. (7)], so that f(0) = 0 and
  f(R) = -g(R) [MIT bag, P(R) = Q(R); g = -Q here] become the natural
  boundary conditions. Only the first spline is non-zero at r = 0
  (B_1(0) = 1), and f(R), g(R) are those of the spline spinors at rmax (see
  add_rmaxBoundary):
  \f[
    S' + \frac{c}{2}[fg]_0^R =
    \frac{c}{2} f(0)^2 \ [\kappa<0] \ {\rm or}\ c^2 f(0)^2 \ [\kappa>0]
    + \frac{c}{4}\left[f(R)^2 - g(R)^2\right] + \frac{c}{2} f(R) g(R),
  \f]
  (the f(0)g(0) terms cancel). The matrix gets its second derivatives. The
  cancellation is exact analytically; on the radial grid it holds only to the
  accuracy of the quadrature at the cavity ends.
  W. R. Johnson, S. A. Blundell, J. Sapirstein, Phys. Rev. A 37, 307 (1988).

  @param Aij     Hamiltonian matrix, 2N x 2N (N f-like splines, then N
                 g-like); updated in place
  @param kappa   Dirac quantum number
  @param alpha   Fine-structure constant
  @param f_rmax  f(rmax) of each spline spinor (see SplineSpinors)
  @param g_rmax  g(rmax) of each spline spinor
*/
void add_JohnsonBoundary(LinAlg::Matrix<double> *Aij, const int kappa,
                         const double alpha, const std::vector<double> &f_rmax,
                         const std::vector<double> &g_rmax);

/*!
  @brief Applies the rmax (MIT bag) boundary term only: f(R) = -g(R).

  @details
  The r = R part of add_JohnsonBoundary (no r = 0 term); also used alone for
  the Hybrid type. With \f$f(R) = \sum_i p_i u_i\f$ and
  \f$g(R) = \sum_i p_i w_i\f$, where \f$u_i\f$, \f$w_i\f$ are f and g of
  spline spinor i at R, adds the second derivatives of
  \f[
    \frac{c}{4}\left[f(R)^2 - g(R)^2\right] + \frac{c}{2} f(R) g(R):
    \quad
    \frac{c}{2}\left(u_iu_j - w_iw_j + u_iw_j + w_iu_j\right).
  \f]
  Without kinetic balance (Johnson), only the last f-like (u = 1) and last
  g-like (w = 1) spinors are non-zero at R. With kinetic balance (Hybrid),
  the kinetic-balance parts of the last two splines of each set also
  contribute, so these are included.

  @param Aij     Hamiltonian matrix (f-like spinors, then g-like); updated in
                 place
  @param alpha   Fine-structure constant
  @param f_rmax  f(rmax) of each spline spinor (see SplineSpinors)
  @param g_rmax  g(rmax) of each spline spinor
*/
void add_rmaxBoundary(LinAlg::Matrix<double> *Aij, const double alpha,
                      const std::vector<double> &f_rmax,
                      const std::vector<double> &g_rmax);

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

  @return Sum minus exact value, for each kappa included. Kappas not present
          in @p basis (e.g., p states for "20s5d") are skipped.
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
