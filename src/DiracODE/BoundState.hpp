#pragma once
#include "AdamsMoulton.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <utility>
#include <vector>
class Grid;

//! Functions and classes used to solve the Dirac equation.
namespace DiracODE {

//==============================================================================
/*!
  @brief Solves bound-state problem for local potential (en < 0).
  @details
  Solves \f$ (H_0 + v - \epsilon_a)F_a = 0 \f$ for the bound state.
  en0 is the initial energy guess (must be reasonably good).
  eps is the convergence target for the energy.
  - @p v is the local potential (e.g., v = v_dir + v_nuc)
  - @p H_off_diag is an optional off-diagonal potential
  - @p alpha: \f$ \alpha = \lambda\alpha_0 \f$ is the effective fine-structure constant
*/
void boundState(DiracSpinor &Fa, const double en0, const std::vector<double> &v,
                const std::vector<double> &H_off_diag = {},
                const double alpha = PhysConst::alpha, double eps = 1.0e-14,
                const DiracSpinor *const VxFa = nullptr,
                const DiracSpinor *const Fa0 = nullptr, double zion = 1,
                double mass = 1.0);

//! Factory overload: constructs a DiracSpinor(n, kappa, gr) and calls boundState().
inline DiracSpinor boundState(
  int n, int kappa, const double en0, const std::shared_ptr<const Grid> &gr,
  const std::vector<double> &v, const std::vector<double> &H_off_diag = {},
  const double alpha = PhysConst::alpha, double eps = 1.0e-14,
  const DiracSpinor *const VxFa = nullptr,
  const DiracSpinor *const Fa0 = nullptr, double zion = 1, double mass = 1.0) {
  DiracSpinor Fnk = DiracSpinor(n, kappa, gr);
  boundState(Fnk, en0, v, H_off_diag, alpha, eps, VxFa, Fa0, zion, mass);
  return Fnk;
}

//! For given energy en, solves DE with correct boundary conditions at the origin
void regularAtOrigin(DiracSpinor &Fa, const double en,
                     const std::vector<double> &v,
                     const std::vector<double> &H_off_diag, const double alpha,
                     const DiracSpinor *const VxFa = nullptr,
                     const DiracSpinor *const Fa0 = nullptr, double zion = 1,
                     double mass = 1.0);

//! For given energy en, solves (local) DE with correct boundary conditions at infinity
void regularAtInfinity(DiracSpinor &Fa, const double en,
                       const std::vector<double> &v,
                       const std::vector<double> &H_off_diag,
                       const double alpha,
                       const DiracSpinor *const VxFa = nullptr,
                       const DiracSpinor *const Fa0 = nullptr, double zion = 1,
                       double mass = 1.0);

namespace Internal {

//==============================================================================
//! Parameters used for Adams-Moulton bound-state solver.
namespace Param {

//! K (# steps) for Adams-Moulton method (between 1 and 12).
constexpr std::size_t K_Adams = 7;
//! Parameter to determine 'asymptotically large r'.
constexpr double cALR = 550.0;
//! Max # attempts at converging bound-state energy.
constexpr int max_its = 99;
//! Fractional size of 'large' energy update steps (~12%).
constexpr double lfrac_de = 0.12;
//! Number of grid points either side of the classical turning point.
constexpr int d_ctp = 4;

//! Order of coefficients in the large-r asymptotic expansion.
constexpr int nx = 15;
//! Convergence threshold for the asymptotic expansion.
constexpr double nx_eps = 1.e-12;

//! Relative cutoff for the large-r tail: extend outwards until |f| drops below
//! tail_cut * max|f|. See @ref extendTail.
constexpr double tail_cut = 1.0e-16;

//! Weighting function for meshing inward/outward solutions at the turning point.
//! Must be positive; index i may be negative [ctp - d_ctp].
constexpr auto weight = [](std::size_t i) {
  return 1.0 / static_cast<double>(i * i + 1);
};

static_assert(
  Param::K_Adams >= 1 && Param::K_Adams <= AdamsMoulton::K_max,
  "\nFAIL in DiracODE: parameter K_Adams must be between 5 and 8\n");

} // namespace Param

//==============================================================================
/*!
  @brief Derivative matrix for the radial Dirac equation, dF/du = D(u)*F(u) + S(u).
  @details
  Implements AdamsMoulton::DerivativeMatrix<std::size_t, double>, using the
  grid index i as the argument type (T = std::size_t). The ODE is solved in
  terms of the grid parameter u (where r = r(u)), so all matrix elements
  include a dr/du Jacobian factor.

  The radial Dirac equation for a central potential v(r) is:

  \f[
    \frac{d}{du}\begin{pmatrix}f\\g\end{pmatrix}
    = \frac{dr}{du}
      \begin{pmatrix}
        -\kappa/r + \alpha H_\text{mag} & \alpha(\varepsilon - v) + 2mc \\
        \alpha(v - \varepsilon)          & \kappa/r - \alpha H_\text{mag}
      \end{pmatrix}
      \begin{pmatrix}f\\g\end{pmatrix}
    + S
  \f]

  where \f$ c = 1/\alpha \f$ is the speed of light in atomic units.
  The optional inhomogeneous source S encodes the exchange interaction:

  \f[
    S_f = -\alpha \, [V_x F_a]_g \, \frac{dr}{du}, \quad
    S_g = +\alpha \, [V_x F_a]_f \, \frac{dr}{du}.
  \f]

  @note Non-copyable; the constructor stores raw pointers to the grid and
        potential arrays, which must outlive this object.
*/
struct DiracDerivative : AdamsMoulton::DerivativeMatrix<std::size_t, double> {

  /*!
    @brief Constructs the Dirac derivative matrix for a given orbital and potential.
    @param in_grid     Radial grid.
    @param in_v        Local potential v(r).
    @param in_k        Orbital kappa quantum number.
    @param in_en       Orbital energy.
    @param in_alpha    Fine-structure constant (alpha).
    @param V_off_diag  Optional off-diagonal (magnetic) potential H_mag(r).
                       If empty, treated as zero.
    @param VxFa        Optional exchange potential VxFa. If nullptr, ignored.
    @param iFa0        Optional inhomogeneous source spinor. If nullptr, ignored.
    @param zion        Effective ionic charge (used for boundary conditions).
    @param in_mass     Effective particle mass in atomic units (default 1).
  */
  DiracDerivative(const Grid &in_grid, const std::vector<double> &in_v,
                  const int in_k, const double in_en, const double in_alpha,
                  const std::vector<double> &V_off_diag = {},
                  const DiracSpinor *const VxFa = nullptr,
                  const DiracSpinor *const iFa0 = nullptr, double zion = 1,
                  double in_mass = 1.0);

  const Grid *const pgr;
  const std::vector<double> *const v;
  const std::vector<double> *const Hmag;
  const DiracSpinor *const VxFa;
  const DiracSpinor *const Fa0;
  const double zion = 1.0;
  const int k;
  const double en, alpha, cc;
  double mass;

  //! D matrix elements (see @ref DiracDerivative for definitions); index i is grid point.
  double a(std::size_t i) const final;
  double b(std::size_t i) const final;
  double c(std::size_t i) const final;
  double d(std::size_t i) const final;
  //! Inhomogeneous source terms from exchange potential VxFa.
  double Sf(std::size_t i) const final;
  double Sg(std::size_t i) const final;

  DiracDerivative(const DiracDerivative &) = delete;
  void operator=(const DiracDerivative &) = delete;
};

//==============================================================================
// Holds the state of the bound-state iteration: the energy brackets used by the
// large (bisection) energy steps, plus the current solution state that is
// passed between the solver steps (trial solution -> energy adjust -> tail).
struct TrackEnGuess {
  // Energy brackets (updated by largeEnergyChange / smallEnergyChangePT):
  int count_toomany = 0;
  int count_toofew = 0;
  double high_en = 0.0;
  double low_en = 0.0;

  // Current solution state, carried between steps:
  double en = 0.0;            // current energy guess
  double eps = 1.0;           // convergence parameter, |de/e|
  std::size_t pinf = 0;       // practical infinity (index after last nonzero)
  std::size_t ctp = 0;        // matching point (near classical turning point)
  int its = 1;                // iteration counter
  bool correct_nodes = false; // current solution has the required # of nodes
  std::vector<double> dg{};   // g_out - g_in near ctp (for the PT energy step)
};

//==============================================================================
//! Builds the trial solution at the current energy: finds pinf and the matching
//! point, integrates out/in, and joins. Returns sofar with pinf, ctp, dg set.
TrackEnGuess trialSolution(DiracSpinor *Fn, const std::vector<double> &v,
                           const std::vector<double> &H_off_diag, double alpha,
                           const DiracSpinor *VxFa, const DiracSpinor *Fa0,
                           double zion, double mass, std::vector<double> &f_in,
                           std::vector<double> &g_in, TrackEnGuess sofar);

//! Counts nodes and updates the energy guess (PT step if # nodes correct, else a
//! large bracketing step). Returns sofar with en, eps, correct_nodes updated.
TrackEnGuess adjustEnergy(DiracSpinor *Fn, int required_nodes, double alpha,
                          TrackEnGuess sofar);

//! Extends the decaying tail past pinf (or, on failure, just zeros the stale
//! tail) and sets Fn.max_pt(). Mutates Fn; reads sofar.
void adjustTail(DiracSpinor *Fn, const std::vector<double> &v,
                const std::vector<double> &H_off_diag, double alpha,
                double zion, double mass, const TrackEnGuess &sofar);

//==============================================================================
//! Returns grid index of "practical infinity": the point where f(r) drops effectively to zero.
std::size_t findPracticalInfinity(const double en, const std::vector<double> &v,
                                  const std::vector<double> &r,
                                  const double alr, const double mass = 1.0);

//! Returns grid index of the classical turning point, where |V(r)| = |E|.
std::size_t findClassicalTurningPoint(const double en,
                                      const std::vector<double> &v,
                                      std::size_t pinf, std::size_t d_ctp);

/*!
  @brief Constructs a trial Dirac solution with correct boundary conditions at both ends.
  @details
  Integrates outwards from the origin and inwards from practical infinity, then
  joins the two solutions at the classical turning point (ctp). If the energy is
  not an eigenvalue, there will be a discontinuity ('kink') in g at ctp. The
  difference dg = g_out - g_in at the joining region is stored and used for the
  perturbation-theory energy update.
*/
void trialDiracSolution(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, std::vector<double> &f_in,
                        std::vector<double> &g_in, const double en,
                        const int ka, const std::vector<double> &v,
                        const std::vector<double> &H_off_diag, const Grid &gr,
                        std::size_t ctp, std::size_t d_ctp, std::size_t pinf,
                        const double alpha,
                        const DiracSpinor *const VxFa = nullptr,
                        const DiracSpinor *const Fa0 = nullptr, double zion = 1,
                        double mass = 1.0);

//! Returns the number of nodes (sign changes) in f up to index maxi.
int countNodes(const std::vector<double> &f, const std::size_t maxi);

//! Makes a large (bisection-style) energy update; updates TrackEnGuess accordingly.
void largeEnergyChange(double *en, TrackEnGuess *sofar, double frac_de,
                       bool toomany_nodes);

//! Returns an updated energy using first-order perturbation theory, given f and dg at ctp.
//! Also tightens the energy brackets in sofar (sign of the PT step tells which
//! side of the eigenvalue en lies on).
double smallEnergyChangePT(const double en, const double anorm,
                           const std::vector<double> &f,
                           const std::vector<double> &dg, std::size_t ctp,
                           std::size_t d_ctp, const double alpha,
                           TrackEnGuess *sofar);

/*!
  @brief Integrates the Dirac equation outwards from the origin.
  @details
  Integrates up to index @p final (not inclusive); if final=0, integrates to
  f.size(). The solution satisfies the boundary condition at r=0 but not at large r.
*/
void solve_Dirac_outwards(std::vector<double> &f, std::vector<double> &g,
                          const DiracDerivative &Hd, std::size_t final = 0);

/*!
  @brief Integrates the Dirac equation inwards from pinf to ctp.
  @details
  The solution satisfies the boundary condition at large r (practical infinity)
  but not at the origin.
*/
void solve_Dirac_inwards(std::vector<double> &f, std::vector<double> &g,
                         const DiracDerivative &Hd, std::size_t ctp,
                         std::size_t pinf, double mass = 1.0);

//==============================================================================
/*!
  @brief Adjusts the practical infinity so the orbital tail ends at the cutoff.
  @details
  After the bound state has been solved (energy converged, in/out solutions
  joined), the orbital is defined up to the "practical infinity" @p pinf, which
  is set by a fixed criterion [findPracticalInfinity] and so does not in general
  sit at a controlled value of |f|. This routine moves @p pinf so that the last
  nonzero point has |f| ~ @p tail_cut * max|f|, i.e. so the orbital ends where it
  has properly decayed (no abrupt truncation, and no redundant near-zero tail).

  - If |f| at @p pinf is still above the cutoff, the decaying tail is continued
    outward (towards larger r). Two methods are used, in order of preference:
    - Preferred: continue the outward Adams-Moulton integration of the
      homogeneous radial Dirac equation, seeded from the last @c K_Adams solved
      points before @p pinf (the same solver used for the bulk of the orbital).
    - Fallback: in the classically-forbidden region the outward solution is
      numerically unstable, as the exponentially growing solution eventually
      dominates. This shows up as f(r) either changing sign or growing in
      magnitude. If that happens, the remaining tail is instead built from the
      large-r asymptotic (exponential-decay) expansion [AsymptoticSpinor],
      scaled to match the solved orbital at @p pinf. This is the same expansion
      used for the inward boundary condition; see the methods document,
      Sec. 3.3.1 (following Johnson 2007).

  - If |f| at @p pinf is already below the cutoff, @p pinf is instead reduced
    (moved inward) to the first point where |f| >= cutoff, and the trimmed tail
    is zeroed.

  Must be called on the converged (but not-yet-normalised) orbital: the cutoff
  is relative to max|f|, which requires the full solution. The orbital is left
  un-normalised; the caller should (re)compute the norm afterwards.

  @param f         Upper radial component, modified in the tail region.
  @param g         Lower radial component, modified in the tail region.
  @param Hd        Homogeneous Dirac derivative matrix (built with the converged
                   energy; exchange/source terms should not be included).
  @param pinf      Current practical infinity (index after final nonzero point).
  @param tail_cut  Relative cutoff: target |f(new_pinf-1)| ~ tail_cut * max|f|.

  @return The new practical infinity (index after the final nonzero point).

  @note Components beyond the returned index are set to zero.
  @warning @p Hd must use the same converged energy as the solved orbital;
           passing an inconsistent energy will produce a discontinuous tail.
*/
std::size_t extendTail(std::vector<double> &f, std::vector<double> &g,
                       const DiracDerivative &Hd, std::size_t pinf,
                       double tail_cut);

/*!
  @brief Joins the inward and outward solutions into a single wavefunction.
  @details
  Produces a solution with correct boundary conditions at both r=0 and large r,
  matched at ctp using a weighted average over [ctp-d_ctp, ctp+d_ctp].
  The resulting wavefunction may not be smooth at the joining point if the
  energy is not an eigenvalue; the discontinuity in g is stored in dg.
*/
void joinInOutSolutions(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg,
                        const std::vector<double> &f_in,
                        const std::vector<double> &g_in, std::size_t ctp,
                        std::size_t d_ctp, std::size_t pinf);

} // namespace Internal
} // namespace DiracODE
