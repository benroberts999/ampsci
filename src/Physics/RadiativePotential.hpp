#pragma once
#include <vector>

//! @brief Ginges-Flambaum QED Radiative Potential
//! @details
//!  - [1] V. V. Flambaum and J. S. M. Ginges, Phys. Rev. A 72, 052115 (2005).
//!  - Including small finite-nuclear size corrections:
//!  - [2] J. S. M. Ginges and J. C. Berengut, Phys. Rev. A 93, 052509 (2016).
//!
//! Functions calculate V_rad(r) for each r, given r_N (nuclear radius) in
//! atomic units, z, and alpha.
//! Note sign: H -> H - V_rad.
namespace RadiativePotential {

//******************************************************************************
//! Small class to hold the radiative potential (for each l)
class Vrad {
public:
  //! Gets Vel = V_euh + V_Se_l + ....  For given l. If l>max, returns V[max]
  const std::vector<double> &get_Vel(const int l = 0) const {
    if (l >= (int)Vel.size())
      return Vel.back();
    return Vel[std::size_t(l)];
  }
  //! Gets Hmag for given l. If l>max, returns H[max]
  const std::vector<double> &get_Hmag(const int l = 0) const {
    if (l >= (int)Hmag.size())
      return Hmag.back();
    return Hmag[std::size_t(l)];
  }

  bool ok() const { return !Vel.empty() && !Hmag.empty() && num_points > 0; }

  void set_Vel(const std::vector<double> &in_v, const int l = 0) {
    ensure_sized(l, in_v.size());
    Vel[std::size_t(l)] = in_v;
  }

  void set_Hmag(const std::vector<double> &in_v, const int l = 0) {
    ensure_sized(l, in_v.size());
    Hmag[std::size_t(l)] = in_v;
  }

private:
  std::size_t num_points = 0;
  std::vector<std::vector<double>> Vel = {{}}; // Vel[l][r], each l
  std::vector<std::vector<double>> Hmag = {{}};

  void ensure_sized(int l, std::size_t in_size) {
    if (l >= (int)Vel.size())
      Vel.resize(std::size_t(l + 1));
    if (l >= (int)Hmag.size())
      Hmag.resize(std::size_t(l + 1));
    if (in_size > num_points) {
      num_points = in_size;
      for (auto &v : Vel)
        v.resize(num_points);
      for (auto &v : Hmag)
        v.resize(num_points);
    }
  }
};

//******************************************************************************
struct RadPot_params {
  double r, rN, z, alpha;
};

// Fitting params From: [Flambaum/Ginges (2005), Ginges/Berengut (2016)]
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
//! Simple V(r) = -(Z^2*alpha) * Exp(-r/lc), lc = hbar / (mc) = alpha
//! @details Note: Does NOT include Ak factor - must be included from input!
//! Does not include FNS effects (rN ignored)
double vSimpleExp(double r, double, double z, double alpha);

//******************************************************************************
//! Uehling (vaccuum polarisation)
double vUehling(double r, double rN, double z, double alpha);

namespace Helper {
// Uehling helper functions:
double vUehcommon(double t, double chi);
double vUehf_smallr(double r, double rN, double chi);
double vUehf_larger(double r, double rN, double chi);
double gslfunc_Ueh_smallr(double t, void *p);
double gslfunc_Ueh_larger(double t, void *p);
} // namespace Helper

//******************************************************************************
//! Self-energy (electric), high-frequency
double vSEh(double r, double rN, double z, double alpha);
//! Self-energy (electric), low-frequency
double vSEl(double r, double rN, double z, double alpha);

namespace Helper {
// Self-energy (elec) helper functions:
double gb_GSEh_smallr(double r, double rN, double chi);
double gb_GSEh_larger(double r, double rN, double chi);
double gb_I1(double t, double z, double alpha);
double gb_I2(double t, double r, double rN, double z, double alpha);
double gslfunc_SEh_smallr(double t, void *p);
double gslfunc_SEh_larger(double t, void *p);
} // namespace Helper

//******************************************************************************
//! Self-energy, magnetic form factor (off-diagonal H function)
double vSE_Hmag(double r, double rN, double z, double alpha);

namespace Helper {
// Self-energy (mag) helper functions:
double vSE_Hmag(double r, double rN, double z, double alpha);
double seJmag(double x, double y);
double gslfunc_SEmag(double t, void *p);
} // namespace Helper

} // namespace RadiativePotential
