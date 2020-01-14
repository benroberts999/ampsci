#pragma once

namespace RadiativePotential {

//******************************************************************************
struct RadPot_params {
  double r, rN, z, alpha;
};

struct Fit_AB {
  // [1] J. S. M. Ginges and J. C. Berengut, Phys. Rev. A 93, 052509 (2016).
  // [1] V. V. Flambaum and J. S. M. Ginges, Phys. Rev. A 72, 052115 (2005).
  double a0 = 1.071, a1 = 0.0, a2 = -1.976, a3 = -2.128, a4 = 0.169;
  double bsp0 = 0.074, bsp1 = 0.35, bsp2 = 0.0;
  double bd0 = 0.056, bd1 = 0.050, bd2 = 0.195;
  //
  double Al(int l, double za) {
    auto x = za - 80.0 / 137.036;
    if (l < 2)
      return a0 + (a1 * x) + (a2 * x * x) + (a3 * x * x * x) +
             (a4 * x * x * x * x);
    return 0.0;
  }
  double Bl(int l, double za) {
    if (l < 2)
      return bsp0 + (bsp1 * za) + (bsp2 * za * za);
    if (l == 2)
      return bd0 + (bd1 * za) + (bd2 * za * za);
    return 0.0;
  }
};

//******************************************************************************
// Uehling (vaccuum pol.)
double vUehling(double r, double rN, double z, double alpha);

// Uehling helper functions:
double vUehcommon(double t, double chi);
double vUehf_smallr(double r, double rN, double chi);
double vUehf_larger(double r, double rN, double chi);
double gslfunc_Ueh_smallr(double t, void *p);
double gslfunc_Ueh_larger(double t, void *p);

//******************************************************************************
// Self-energy (electric)
double vSEh(double r, double rN, double z, double alpha);
double vSEl(double r, double rN, double z, double alpha);

// Self-energy helper functions:
double gb_GSEh_smallr(double r, double rN, double chi);
double gb_GSEh_larger(double r, double rN, double chi);
double gb_I1(double t, double z, double alpha);
double gb_I2(double t, double r, double rN, double z, double alpha);
double gslfunc_SEh_smallr(double t, void *p);

} // namespace RadiativePotential
