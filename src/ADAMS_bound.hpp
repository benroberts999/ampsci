#pragma once
#include <array>
#include <vector>
class DiracSpinor;
class Grid;

namespace ADAMS {

const int AMO = 8;
static_assert(
    AMO >= 5 && AMO <= 8,
    "\nFAIL 8 in Adams (.h): parameter AMO must be between 5 and 8\n");

void solveDBS(DiracSpinor &psi, const std::vector<double> &v, const Grid &rgrid,
              const double alpha, int log_dele = 0);

int findPracticalInfinity(const double en, const std::vector<double> &v,
                          const std::vector<double> &r, const double alr);

int findClassicalTurningPoint(const double en, const std::vector<double> &v,
                              const int pinf, const int d_ctp);

void trialDiracSolution(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, const double en, const int ka,
                        const std::vector<double> &v, const Grid &gr,
                        const int ctp, const int d_ctp, const int pinf,
                        const double alpha);

int countNodes(const std::vector<double> &f, const int maxi);

void largeEnergyChange(double &en, int &count_toomany, int &count_toofew,
                       double &high_en, double &low_en, double lfrac_de,
                       bool toomany_nodes);

double calcNorm(const std::vector<double> &f, const std::vector<double> &g,
                const std::vector<double> &drdt, const double h,
                const int pinf);

double smallEnergyChangePT(const double en, const double anorm,
                           const std::vector<double> &f,
                           const std::vector<double> &dg, int ctp, int d_ctp,
                           double alpha, int count_toofew, int count_toomany,
                           double low_en, double high_en);

void outwardAM(std::vector<double> &p, std::vector<double> &q, const double en,
               const std::vector<double> &v, const int ka,
               const std::vector<double> &r, const std::vector<double> &drdt,
               const double h, const int ctp, const double alpha);

void inwardAM(std::vector<double> &f, std::vector<double> &g, const double en,
              const std::vector<double> &v, const int ka,
              const std::vector<double> &r, const std::vector<double> &drdt,
              const double h, const int ctp, const int pinf,
              const double alpha);

void adamsMoulton(std::vector<double> &f, std::vector<double> &g,
                  const double en, const std::vector<double> &v, const int ka,
                  const std::vector<double> &r, const std::vector<double> &drdu,
                  const double du, const int ni, const int nf,
                  const double alpha);

void joinInOutSolutions(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg,
                        const std::vector<double> &f_in,
                        const std::vector<double> &g_in, const int ctp,
                        const int d_ctp, const int pinf);

//******************************************************************************
template <int N> struct AdamsCoefs {};
//------------------------------------------------------------------------------
template <> struct AdamsCoefs<8> {
  // Adams coefs:
  static const int N = 8;
  const int AMa[N] = {-33953,   312874,  -1291214, 3146338,
                      -5033120, 5595358, -4604594, 4467094};
  const double AMd = 1. / 3628800;
  const int AMaa = 1070017.;
  // Outward coefs:
  const int OIe[N][N] = {
      {-1338, 2940, -2940, 2450, -1470, 588, -140, 15},
      {-240, -798, 1680, -1050, 560, -210, 48, -5},
      {60, -420, -378, 1050, -420, 140, -30, 3},
      {-32, 168, -672, 0, 672, -168, 32, -3},
      {30, -140, 420, -1050, 378, 420, -60, 5},
      {-48, 210, -560, 1050, -1680, 798, 240, -15},
      {140, -588, 1470, -2450, 2940, -2940, 1338, 105},
      {-960, 3920, -9408, 14700, -15680, 11760, -6720, 2283}};
  const int OIa[N] = {-105, 15, -5, 3, -3, 5, -15, 105};
  const int OId = 840.;
};
//------------------------------------------------------------------------------
template <> struct AdamsCoefs<7> {
  // Adams coefs:
  static const int N = 7;
  const int AMa[N] = {1375, -11351, 41499, -88547, 123133, -121797, 139849};
  const double AMd = 1. / 120960;
  const int AMaa = 36799;
  // Outward coefs:
  const int OIe[N][N] = {{-609, 1260, -1050, 700, -315, 84, -10},
                         {-140, -329, 700, -350, 140, -35, 4},
                         {42, -252, -105, 420, -126, 28, -3},
                         {-28, 126, -420, 105, 252, -42, 4},
                         {35, -140, 350, -700, 329, 140, -10},
                         {-84, 315, -700, 1050, -1260, 609, 60},
                         {490, -1764, 3675, -4900, 4410, -2940, 1089}};
  const int OIa[N] = {-60, 10, -4, 3, -4, 10, -60};
  const int OId = 420;
};
//------------------------------------------------------------------------------
template <> struct AdamsCoefs<6> {
  // Adams coefs:
  static const int N = 6;
  const int AMa[N] = {-863, 6312, -20211, 37504, -46461, 65112};
  const double AMd = 1. / 60480;
  const int AMaa = 19087;
  // Outward coefs:
  const int OIe[N][N] = {
      {-77, 150, -100, 50, -15, 2}, {-24, -35, 80, -30, 8, -1},
      {9, -45, 0, 45, -9, 1},       {-8, 30, -80, 35, 24, -2},
      {15, -50, 100, -150, 77, 10}, {-72, 225, -400, 450, -360, 147}};
  const int OIa[N] = {-10, 2, -1, 1, -2, 10};
  const int OId = 60;
};
//------------------------------------------------------------------------------
template <> struct AdamsCoefs<5> {
  // Adams coefs:
  static const int N = 5;
  const int AMa[N] = {27, -173, 482, -798, 1427};
  const double AMd = 1. / 1440;
  const int AMaa = 475;
  // Outward coefs:
  const int OIe[N][N] = {{-65, 120, -60, 20, -3},
                         {-30, -20, 60, -15, 2},
                         {15, -60, 20, 30, -3},
                         {-20, 60, -120, 65, 12},
                         {75, -200, 300, -300, 137}};
  const int OIa[N] = {-12, 3, -2, 3, -12};
  const int OId = 60;
};

} // namespace ADAMS
