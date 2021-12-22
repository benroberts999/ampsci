#include "Physics/NuclearPotentials.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "NuclearData.hpp" //for Isotope
#include "Physics/AtomData.hpp"
#include "Physics/NuclearData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "qip/String.hpp"
#include <cmath>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <string>
#include <vector>

namespace Nuclear {

//******************************************************************************
Type parseType(const std::string &str_type) {
  if (qip::ci_wildcard_compare(str_type, "fermi"))
    return Type::Fermi;
  if (qip::ci_wildcard_compare(str_type, "spher*") ||
      qip::ci_wildcard_compare(str_type, "ball"))
    return Type::spherical;
  if (qip::ci_wildcard_compare(str_type, "point*"))
    return Type::point;
  if (qip::ci_wildcard_compare(str_type, "gaus*"))
    return Type::Gaussian;
  std::cout << "\n⚠️  WARNING: Unkown nucleus type: " << str_type
            << "; defaulting to Fermi\n"
            << "Options are: Fermi, spherical, pointlike, Gaussian\n\n";
  return Type::Fermi;
}
std::string parseType(Type type) {
  if (type == Type::Fermi)
    return "Fermi";
  if (type == Type::spherical)
    return "spherical";
  if (type == Type::point)
    return "pointlike";
  if (type == Type::Gaussian)
    return "Gaussian";
  return "Fermi";
}

//******************************************************************************
Parameters::Parameters(int in_z, int in_a, const std::string &str_type,
                       double in_rrms, double in_t)
    : z(in_z),
      a((in_a < 0) ? AtomData::defaultA(z) : in_a),
      type(parseType(str_type)),
      r_rms(in_rrms),
      t(in_t <= 0 ? approximate_t_skin(a) : in_t) {
  if (r_rms < 0) {
    r_rms = find_rrms(z, a);
    if (r_rms <= 0)
      r_rms = approximate_r_rms(a);
  }
}
//------------
Parameters::Parameters(const std::string &z_str, int in_a,
                       const std::string &str_type, double in_rrms, double in_t)
    : z(AtomData::atomic_Z(z_str)),
      a((in_a < 0) ? AtomData::defaultA(z) : in_a),
      type(parseType(str_type)),
      r_rms(in_rrms),
      t(in_t <= 0 ? approximate_t_skin(a) : in_t) {
  if (r_rms < 0) {
    r_rms = find_rrms(z, a);
    if (r_rms <= 0)
      r_rms = approximate_r_rms(a);
  }
}

//******************************************************************************
std::vector<double> sphericalNuclearPotential(double Z, double rnuc,
                                              const std::vector<double> &rgrid)
// Potential due to a spherical nucleus, with (charge) radius, rnuc.
// Note: rnuc must be given in "fermi" (fm, femto-metres).
// rnuc = 0 corresponds to zeroNucleus
{
  std::vector<double> vnuc;
  vnuc.reserve(rgrid.size());

  // Fill the vnuc array with spherical nuclear potantial
  const double rN = rnuc / PhysConst::aB_fm; // convert fm -> au
  const double rn2 = rN * rN;
  const double rn3 = rn2 * rN;
  for (auto r : rgrid) {
    const double temp_v =
        (r < rN) ? Z * (r * r - 3.0 * rn2) / (2.0 * rn3) : -Z / r;
    vnuc.push_back(temp_v);
  }

  return vnuc;
}

//******************************************************************************
std::vector<double> GaussianNuclearPotential(double Z, double r_rms,
                                             const std::vector<double> &rgrid)

{
  std::vector<double> vnuc;
  vnuc.reserve(rgrid.size());

  // Fill the vnuc array with Gaussian nuclear potantial
  const double rN = r_rms / PhysConst::aB_fm; // convert fm -> au
  const auto k = std::sqrt(3.0 / 2);
  for (auto r : rgrid) {
    vnuc.push_back(-Z * std::erf(k * r / rN) / r);
  }

  return vnuc;
}

//******************************************************************************
std::vector<double> fermiNuclearPotential(double z, double t, double c,
                                          const std::vector<double> &rgrid)
// Uses a Fermi-Dirac distribution for the nuclear potential.
//
// rho(r) = rho_0 {1 + Exf[(r-c)/a]}^-1
// V(r) = -(4 Pi)/r [A+B]
//   A = Int[ rho(x) x^2 , {x,0,r}]
//   B = r * Int[ rho(x) x , {x,r,infty}]
// rho_0 is found by either:
//   * V(infinity) = -Z/r , or equivilantly
//   * int rho(r) d^3r = Z
//
// Depends on:
//   * t: skin thickness [90 to 10% fall-off range]
//     note: t = a[4 ln(3)]
//   * c: half-density radius (see above for link c and r_rms)
//
// t and c are input values. In 'fermi' or fm (fempto metres)
//
// V(r) is expressed in terms of Complete Fermi-Dirac intagrals.
// These are computed using the GSL libraries.
// gnu.org/software/gsl/manual/html_node/Complete-Fermi_002dDirac-Integrals
{
  std::vector<double> vnuc;
  vnuc.reserve(rgrid.size());

  const double a = t / FourLn3;
  const double coa = c / a;
  // Use GSL for the Complete Fermi-Dirac Integrals:
  const double f2 = gsl_sf_fermi_dirac_2(coa);
  for (auto r : rgrid) {
    double t_v = -z / r;
    const double roa = PhysConst::aB_fm * r / a; // convert fm <-> atomic
    const double roc = r / c * PhysConst::aB_fm;
    if (roc < 10.0) {
      const auto coa2 = coa * coa;
      const auto roa3 = roa * roa * roa;
      const double xf1 = gsl_sf_fermi_dirac_1(roa - coa);
      const double xf2 = gsl_sf_fermi_dirac_2(roa - coa);
      const double tX = -roa3 - 2.0 * coa * (Pi2 + coa2) +
                        roa * (Pi2 + 3.0 * coa2) + 6.0 * roa * xf1 - 12.0 * xf2;
      t_v += t_v * tX / (12.0 * f2);
    }
    vnuc.push_back(t_v);
  }
  return vnuc;
}

//******************************************************************************
std::vector<double> fermiNuclearDensity_tcN(double t, double c, double Z_norm,
                                            const Grid &grid)
// Integrate[ rho(r) , dV ] = Integrate[ 4pi * r^2 * rho(r) , dr ] = Z_norm
// Znorm = Z for nuclear chare density; Z_norm = 1 for nuclear density.
{
  std::vector<double> rho;
  rho.reserve(grid.num_points());

  // form un-normalised rho:
  const double a = t / FourLn3;
  const double coa = c / a;
  for (const auto r : grid.r()) {
    const double roa = PhysConst::aB_fm * r / a;
    const double roc = r / c * PhysConst::aB_fm;
    if (roc < 10.0) {
      rho.push_back(1.0 / (1.0 + std::exp(roa - coa)));
    } else {
      rho.push_back(0.0);
    }
  }

  // Find rho0, normalisation constant + re-scale (normalise rho)
  const double volume_integral =
      NumCalc::integrate(grid.du(), 0, 0, grid.r(), grid.r(), rho,
                         grid.drdu()) *
      4.0 * M_PI;
  const double rho0 = Z_norm / volume_integral;
  for (auto &rhoi : rho)
    rhoi *= rho0;

  return rho;
}

//******************************************************************************
std::vector<double> formPotential(const Parameters &params,
                                  const std::vector<double> &r) {
  const auto z = params.z;
  const auto nucleus_type = params.type;
  const auto r_rms = params.r_rms;
  const auto t = params.t;

  switch (nucleus_type) {

  case Nuclear::Type::Fermi: {
    const auto chdr = Nuclear::c_hdr_formula_rrms_t(r_rms, t);
    return Nuclear::fermiNuclearPotential(z, t, chdr, r);
  }
  case Nuclear::Type::spherical: {
    const auto r_n = Nuclear::c_hdr_formula_rrms_t(r_rms, 0.0); // right?
    return Nuclear::sphericalNuclearPotential(z, r_n, r);
  }
  case Nuclear::Type::point:
    return Nuclear::sphericalNuclearPotential(z, 0.0, r);
  case Nuclear::Type::Gaussian:
    return Nuclear::GaussianNuclearPotential(z, r_rms, r);

  default:
    std::cerr << "\nFAIL Nuclear:266 - invalid nucleus type?\n";
    std::abort();
  }
}

} // namespace Nuclear
