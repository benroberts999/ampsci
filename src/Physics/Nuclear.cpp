#include "Nuclear.hpp"  //:(
#include "AtomInfo.hpp" //:(
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Nuclear_DataTable.hpp"
#include "PhysConst_constants.hpp"
#include <cmath>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <string>
#include <vector>

namespace Nuclear {

//******************************************************************************
Type parseType(const std::string &str_type) {
  if (str_type == "Fermi")
    return Type::Fermi;
  if (str_type == "spherical")
    return Type::spherical;
  if (str_type == "zero")
    return Type::point;
  if (str_type == "point")
    return Type::point;
  if (str_type == "pointlike")
    return Type::point;
  return Type::Fermi;
}

//******************************************************************************
Parameters::Parameters(int in_z, int in_a, std::string str_type, double in_rrms,
                       double in_t)
    : z(in_z),                                      //
      a((in_a < 0) ? AtomInfo::defaultA(z) : in_a), //
      type(parseType(str_type)),                    //
      t(in_t <= 0 ? approximate_t_skin(a) : in_t),  //
      r_rms(in_rrms)                                //
{
  if (r_rms < 0) {
    r_rms = find_rrms(z, a);
    if (r_rms <= 0)
      r_rms = approximate_r_rms(a);
  }
}

//******************************************************************************
Isotope findIsotopeData(int z, int a) {
  for (const auto &nucleus : NuclearDataTable) {
    if (nucleus.Z == z && nucleus.A == a)
      return nucleus;
  }
  return Isotope{z, a, -1, 0, 0, -1};
}

std::vector<Isotope> findIsotopeList(int z) {
  std::vector<Isotope> isotopes;
  for (const auto &nucleus : NuclearDataTable) {
    if (nucleus.Z == z)
      isotopes.push_back(nucleus);
  }
  return isotopes;
}

double find_rrms(int z, int a) {
  auto nuc = findIsotopeData(z, a);
  if (!nuc.r_ok())
    return 0;
  return nuc.r_rms;
}

double find_mu(int z, int a) {
  auto nuc = findIsotopeData(z, a);
  if (!nuc.mu_ok())
    return 0;
  return nuc.mu;
}

int find_parity(int z, int a) {
  auto nuc = findIsotopeData(z, a);
  if (!nuc.parity_ok())
    return 0;
  return nuc.parity;
}

double find_spin(int z, int a) {
  auto nuc = findIsotopeData(z, a);
  if (!nuc.I_ok())
    return -1.0;
  return nuc.I_N;
}

//******************************************************************************
double approximate_r_rms(int A)
// Returns approximate root-mean-square charge radius in fm [1.e-15 m]
// https://www.sciencedirect.com/science/article/pii/S0092640X12000265
// https://www-nds.iaea.org/radii/
{
  return (A < 10) ? 1.15 * std::pow(A, 0.333)
                  : 0.836 * std::pow(A, 0.333) + 0.570;
}

//******************************************************************************
double c_hdr_formula_rrms_t(double rrms, double t)
// Calculates half-density radius, given rms charge radius, and t.
// Formula from Ginges, Volotka, Fritzsche, Phys. Rev. A 96, 1 (2017).
// 4 ln(3) = 4.394449155, pi^2 = 9.86960440
{
  const double a = t / 4.394449155;
  if (rrms < t) {
    // this is little dodgy? but formula prob only works large A
    return std::sqrt((5. / 3) * rrms * rrms);
  }
  return std::sqrt((5. / 3) * rrms * rrms - (7. / 3) * (9.86960440 * a * a));
}

//******************************************************************************
double rrms_formula_c_t(double c, double t)
// Calculates  rms charge radius, given half-density radius (c), and t.
// Formula from Ginges, Volotka, Fritzsche, Phys. Rev. A 96, 1 (2017).
// 4 ln(3) = 4.39445, pi^2 = 9.87
{
  const double a = t / 4.394449155;
  return std::sqrt(0.2 * (3.0 * c * c + 7.0 * a * a * 9.86960440));
}

//******************************************************************************

constexpr double approximate_t_skin(int) { return default_t; }

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
    double temp_v = (r < rN) ? Z * (r * r - 3.0 * rn2) / (2.0 * rn3) : -Z / r;
    vnuc.push_back(temp_v);
  }

  return vnuc;
}

//******************************************************************************
std::vector<double> fermiNuclearPotential(double Z, double t, double c,
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
//   * c: half-density raius [rho(c)=0.5 rho0]
//
// t and c are input values. In 'fermi' of fm (fempto metres)
//
// V(r) is expressed in terms of Complete Fermi-Dirac intagrals.
// These are computed using the GSL libraries.
// gnu.org/software/gsl/manual/html_node/Complete-Fermi_002dDirac-Integrals
{
  std::vector<double> vnuc;
  vnuc.reserve(rgrid.size());

  const double a = 0.2275598067 * t; // t = a*[4 ln(3)]
  const double coa = c / a;
  // Use GSL for the Complete Fermi-Dirac Integrals:
  const double F2 = gsl_sf_fermi_dirac_2(coa);
  const double pi2 = M_PI * M_PI;
  for (auto r : rgrid) {
    double t_v = -Z / r;
    const double roa = PhysConst::aB_fm * r / a; // convert fm <-> atomic
    const double roc = r / c * PhysConst::aB_fm;
    if (roc < 10.0) {
      double coa2 = coa * coa;
      double xF1 = gsl_sf_fermi_dirac_1(roa - coa);
      double xF2 = gsl_sf_fermi_dirac_2(roa - coa);
      double tX = -std::pow(roa, 3) - 2 * coa * (pi2 + coa2) +
                  roa * (pi2 + 3 * coa2) + 6 * roa * xF1 - 12 * xF2;
      t_v += t_v * tX / (12.0 * F2);
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
  rho.reserve(grid.ngp);

  const double a = 0.2275598067 * t;
  const double coa = c / a;
  for (auto r : grid.r) {
    const double roa = PhysConst::aB_fm * r / a;
    const double roc = r / c * PhysConst::aB_fm;
    if (roc < 10.0) {
      rho.push_back(1.0 / (1.0 + std::exp(roa - coa)));
    } else {
      rho.push_back(0.0);
    }
  }

  double Norm =
      NumCalc::integrate({&grid.r, &grid.r, &rho, &grid.drdu}, grid.du) * 4.0 *
      M_PI;
  double rho0 = Z_norm / Norm;

  for (auto &rhoi : rho)
    rhoi *= rho0;

  return rho;
}

//******************************************************************************
std::vector<double> formPotential(Parameters params, int z, int,
                                  const std::vector<double> &r) {
  const auto nucleus_type = params.type;
  const auto r_rms = params.r_rms;
  const auto t = params.t;

  switch (nucleus_type) {

  case Nuclear::Type::Fermi: {
    auto chdr = Nuclear::c_hdr_formula_rrms_t(r_rms, t);
    return Nuclear::fermiNuclearPotential(z, t, chdr, r);
  }
  case Nuclear::Type::spherical: {
    auto r_n = Nuclear::c_hdr_formula_rrms_t(r_rms, 0.0); // right?
    return Nuclear::sphericalNuclearPotential(z, r_n, r);
  }
  case Nuclear::Type::point:
    return Nuclear::sphericalNuclearPotential(z, 0.0, r);

  default:
    std::cerr << "\nFAIL Nuclear:266 - invalid nucleus type?\n";
    std::abort();
  }
}

} // namespace Nuclear
