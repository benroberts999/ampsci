#pragma once
#include "AdamsMoulton.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <utility>
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {

//==============================================================================
//! @brief Solves bound-state problem for local potential (en < 0)
/*! @details
\f[ (H_0 + v - \epsilon_a)F_a = 0\f]
en0 is initial energy guess (must be reasonably good).
log_eps: log10(eps); eps is convergence target for energy.
  - v is local potential (e.g., v = v_dir + v_nuc)
  - H_off_diag is optional off-diagonal potential.
  - alpha: \f$\alpha = \lambda\alpha_0\f$ is the effective value of
fine-structure constant
*/
void boundState(DiracSpinor &Fa, const double en0, const std::vector<double> &v,
                const std::vector<double> &H_off_diag = {},
                const double alpha = PhysConst::alpha, double eps = 1.0e-14,
                const DiracSpinor *const VxFa = nullptr,
                const DiracSpinor *const Fa0 = nullptr, double zion = 1);

//! For given energy en, solves DE with correct boundary conditions at the origin
void regularAtOrigin(DiracSpinor &Fa, const double en,
                     const std::vector<double> &v,
                     const std::vector<double> &H_off_diag, const double alpha,
                     const DiracSpinor *const VxFa = nullptr,
                     const DiracSpinor *const Fa0 = nullptr, double zion = 1);

//! For given energy en, solves (local) DE with correct boundary conditions at infinity
void regularAtInfinity(DiracSpinor &Fa, const double en,
                       const std::vector<double> &v,
                       const std::vector<double> &H_off_diag,
                       const double alpha,
                       const DiracSpinor *const VxFa = nullptr,
                       const DiracSpinor *const Fa0 = nullptr, double zion = 1);

namespace Internal {

//==============================================================================
// Parameters used for Adams-Moulton mehtod:
namespace Param {

// K (# steps) for Adams-Moulton method (between 1 and 12)
constexpr std::size_t K_Adams = 7;
// Parameter to determine 'assymptotically large r [kinda..]' (=800)
constexpr double cALR = 550.0;
// Max # attempts at converging [sove bs] (30)
constexpr int max_its = 99;
// 'large' energy variations (0.1 => 10%)
constexpr double lfrac_de = 0.12;
// Num points past ctp +/- d_ctp.
constexpr int d_ctp = 4;

// order of the expansion coeficients in large-r asymptotic expansion  (15 orig.)
constexpr int nx = 15;
// convergance for expansion in asymptotic expansion (10^-8)
constexpr double nx_eps = 1.e-12;

// weighting function for meshing in/out solutions
// nb: must be positive, but i may be negative [ctp - d_ctp]
constexpr auto weight = [](std::size_t i) {
  return 1.0 / static_cast<double>(i * i + 1);
};

static_assert(
    Param::K_Adams >= 1 && Param::K_Adams <= AdamsMoulton::K_max,
    "\nFAIL in DiracODE: parameter K_Adams must be between 5 and 8\n");

} // namespace Param

//==============================================================================
//! Matrix which defines Dirac derivative: (dF/dr) = D*F
struct DiracDerivative : AdamsMoulton::DerivativeMatrix<std::size_t, double> {

  DiracDerivative(const Grid &in_grid, const std::vector<double> &in_v,
                  const int in_k, const double in_en, const double in_alpha,
                  const std::vector<double> &V_off_diag = {},
                  const DiracSpinor *const VxFa = nullptr,
                  const DiracSpinor *const iFa0 = nullptr, double zion = 1);
  const Grid *const pgr;
  const std::vector<double> *const v;
  const std::vector<double> *const Hmag;
  const DiracSpinor *const VxFa;
  const DiracSpinor *const Fa0;
  const double zion = 1.0;
  const int k;
  const double en, alpha, cc;

  double a(std::size_t i) const final;
  double b(std::size_t i) const final;
  double c(std::size_t i) const final;
  double d(std::size_t i) const final;
  double Sf(std::size_t i) const final;
  double Sg(std::size_t i) const final;

  DiracDerivative(const DiracDerivative &) = delete;
  void operator=(const DiracDerivative &) = delete;
};

//==============================================================================
// To keep track of current/previous energy guesses
struct TrackEnGuess {
  // Number of times there was too many/too few nodes:
  int count_toomany = 0;
  int count_toofew = 0;
  // Upper and lower energy window before correct # nodes
  double high_en = 0.0;
  double low_en = 0.0;
};

//==============================================================================
// Finds "practical infinity" index, where f(r) drops effectively to zero
std::size_t findPracticalInfinity(const double en, const std::vector<double> &v,
                                  const std::vector<double> &r,
                                  const double alr);

// Finds "classical turning point" index, where |V(r)| = |E|
std::size_t findClassicalTurningPoint(const double en,
                                      const std::vector<double> &v,
                                      std::size_t pinf, std::size_t d_ctp);

// Finds a trial Dirac solution for given energy that has correct boundary conditions.
/*
If it's not an exact solution, there will be a 'kink' in the wavefunction 
around ctp (classical turning point), where the inward and outward solutons 
were joind. The difference between g for these, dg=(gout-gin), is stored, and 
is used for PT
*/
void trialDiracSolution(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, const double en, const int ka,
                        const std::vector<double> &v,
                        const std::vector<double> &H_off_diag, const Grid &gr,
                        std::size_t ctp, std::size_t d_ctp, std::size_t pinf,
                        const double alpha,
                        const DiracSpinor *const VxFa = nullptr,
                        const DiracSpinor *const Fa0 = nullptr,
                        double zion = 1);

// Counts the nodes a given function f
int countNodes(const std::vector<double> &f, const std::size_t maxi);

// Makes a large update to energy; updates 'TrackEnGuess'
void largeEnergyChange(double *en, TrackEnGuess *sofar, double frac_de,
                       bool toomany_nodes);

// Makes a small update to energy using perturbation theory, given f and dg
double smallEnergyChangePT(const double en, const double anorm,
                           const std::vector<double> &f,
                           const std::vector<double> &dg, std::size_t ctp,
                           std::size_t d_ctp, const double alpha,
                           const TrackEnGuess &sofar);

// Solves Dirac equation by integrating outwards from zero.
// Integrates only to 'final' (not inclusive). If final=0, goes to f.size()
// Solution has correct boundary condition at r=0, but not at large r.
void solve_Dirac_outwards(std::vector<double> &f, std::vector<double> &g,
                          const DiracDerivative &Hd, std::size_t final = 0);

// Solves Dirac equation by integrating inwards from 'pinf' to 'ctp'
// Solution has correct boundary condition at large r, but not at small r.
void solve_Dirac_inwards(std::vector<double> &f, std::vector<double> &g,
                         const DiracDerivative &Hd, std::size_t ctp,
                         std::size_t pinf);

// Meshes the two solutions from inwards/outwards integration.
// Produces solution that has correct boundary conditions at 0 and infinity,
// but may not be smooth at the joining point.
void joinInOutSolutions(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg,
                        const std::vector<double> &f_in,
                        const std::vector<double> &g_in, std::size_t ctp,
                        std::size_t d_ctp, std::size_t pinf);

} // namespace Internal
} // namespace DiracODE
