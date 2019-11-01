#pragma once
#include "Adams_coefs.hpp"
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {
namespace Adams {

//******************************************************************************
// Parameters used for Adams-Moulton mehtod:
namespace Param {
constexpr int AMO = 5;            // Adams-Moulton (between 5 and 8)
constexpr AdamsCoefs<AMO> AMcoef; // Adamns-Moulton coeficients [defined .h]
constexpr double cALR = 875;      // 'assymptotically large r [kinda..]' (=800)
constexpr int max_its = 99;       // Max # attempts at converging [sove bs] (30)
constexpr double lfrac_de = 0.12; // 'large' energy variations (0.1 => 10%)
constexpr int d_ctp = 4;          // Num points past ctp +/- d_ctp.

// order of the expansion coeficients in 'inwardAM'  (15 orig.)
constexpr int nx = 15;
// convergance for expansion in `inwardAM'' (10^-8)
constexpr double nx_eps = 1.e-14;
// # of outdir runs [finds first Param::num_loops*AMO+1 points (3)]
constexpr int num_loops = 1;

// weighting function for meshing in/out solutions
// nb: must be positive, but i may be negative (?) [ctp - d_ctp]
constexpr auto weight = [](const int i) {
  return 1.0 / static_cast<double>(i * i + 1);
};

static_assert(
    Param::AMO >= 5 && Param::AMO <= 8,
    "\nFAIL 8 in Adams (.h): parameter AMO must be between 5 and 8\n");

} // namespace Param

//******************************************************************************
class DiracMatrix {
  // Notation:
  // df = af - bg
  // dg = -cf + dg
public:
  DiracMatrix(const Grid &in_grid, const std::vector<double> &in_v,
              const int in_k, const double in_en, const double in_alpha);

  const Grid *const pgr;
  const std::vector<double> *const v;
  const int k;
  const double en, alpha, c2;

  // update a and d for off-diag additional potential (magnetic form-fac, QED)
  double a(std::size_t i) const;
  double b(std::size_t i) const;
  double c(std::size_t i) const;
  double d(std::size_t i) const;
};

struct TrackEnGuess {
  // Number of times there was too many/too few nodes:
  int count_toomany = 0;
  int count_toofew = 0;
  // Upper and lower energy window before correct # nodes
  double high_en = 0.0;
  double low_en = 0.0;
};

// -----------------------------------------------------------------------------
int findPracticalInfinity(const double en, const std::vector<double> &v,
                          const std::vector<double> &r, const double alr);

int findClassicalTurninum_pointsoint(const double en,
                                     const std::vector<double> &v,
                                     const int pinf, const int d_ctp);

void trialDiracSolution(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, const double en, const int ka,
                        const std::vector<double> &v, const Grid &gr,
                        const int ctp, const int d_ctp, const int pinf,
                        const double alpha);

int countNodes(const std::vector<double> &f, const int maxi);

// c++17: could use structured binding, move in/outs to returns
void largeEnergyChange(double *en, TrackEnGuess *sofar, double frac_de,
                       bool toomany_nodes);

double calcNorm(const std::vector<double> &f, const std::vector<double> &g,
                const std::vector<double> &drdt, const double h,
                const int pinf);

double smallEnergyChangePT(const double en, const double anorm,
                           const std::vector<double> &f,
                           const std::vector<double> &dg, const int ctp,
                           const int d_ctp, const double alpha,
                           const TrackEnGuess &sofar);

void outwardAM(std::vector<double> &f, std::vector<double> &g,
               const DiracMatrix &Hd, const int final);

void inwardAM(std::vector<double> &f, std::vector<double> &g,
              const DiracMatrix &Hd, const int ctp, const int pinf);

void adamsMoulton(std::vector<double> &f, std::vector<double> &g,
                  const DiracMatrix &Hd, const int ni, const int nf);

void joinInOutSolutions(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg,
                        const std::vector<double> &f_in,
                        const std::vector<double> &g_in, const int ctp,
                        const int d_ctp, const int pinf);

} // namespace Adams
} // namespace DiracODE
