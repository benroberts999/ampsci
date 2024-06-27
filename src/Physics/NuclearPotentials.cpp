#include "Physics/NuclearPotentials.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "NuclearData.hpp" //for Isotope
#include "Physics/AtomData.hpp"
#include "Physics/NuclearData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "qip/String.hpp"
#include "qip/Vector.hpp"
#include <cmath>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <iostream>
#include <string>
#include <vector>

namespace Nuclear {

//==============================================================================
ChargeDistro parseType(const std::string &str_type) {
  if (qip::ci_wc_compare(str_type, "fermi"))
    return ChargeDistro::Fermi;
  if (qip::ci_wc_compare(str_type, "spher*") ||
      qip::ci_wc_compare(str_type, "ball") ||
      qip::ci_wc_compare(str_type, "uniform"))
    return ChargeDistro::spherical;
  if (qip::ci_wc_compare(str_type, "point*") ||
      qip::ci_wc_compare(str_type, "zero*"))
    return ChargeDistro::point;
  if (qip::ci_wc_compare(str_type, "gaus*"))
    return ChargeDistro::Gaussian;
  if (qip::ci_wc_compare(str_type, "custom"))
    return ChargeDistro::custom;
  std::cout << "\n!!  WARNING: Unknown nucleus type: " << str_type << "\n";

  // spelling suggestions:
  std::cout << "Did you mean: "
            << *qip::ci_closest_match(str_type, {"Fermi", "spherical", "ball",
                                                 "pointlike", "Gaussian", "custom"})
            << " ?\n ";
  std::cout << "Options are: Fermi, spherical (or ball, uniform), pointlike, "
               "Gaussian, or custom (requires readable input file)\n\n";

  return ChargeDistro::Error;
}

//------------------------------------------------------------------------------
std::string parseType(ChargeDistro type) {
  if (type == ChargeDistro::Fermi)
    return "Fermi";
  if (type == ChargeDistro::spherical)
    return "spherical";
  if (type == ChargeDistro::point)
    return "pointlike";
  if (type == ChargeDistro::Gaussian)
    return "Gaussian";
  if (type == ChargeDistro::custom)
    return "custom";
  return "Error";
}

//==============================================================================
Nucleus::Nucleus(int tz, int ta, const std::string &str_type, double trrms,
                 double tt, const std::vector<double> &in_params,
                 const std::string & custom_pot_file_name)
    : m_iso(findIsotopeData(tz, ta < 0 ? AtomData::defaultA(tz) : ta)),
      m_type(parseType(str_type)),
      m_t(tt),
      m_beta(0.0),
      m_custom_pot_file_name(custom_pot_file_name),
      m_params(in_params) {

  if (m_type != ChargeDistro::Fermi)
    m_t = 0.0;

  if (trrms >= 0.0)
    set_rrms(trrms);

  if (ta == 0)
    set_rrms(0.0);

  if (r_rms() <= 0.0 && trrms < 0.0 && ta != 0) {
    const auto approx_rrms = approximate_r_rms(a(), z());
    set_rrms(approx_rrms);
    std::cout << "\n\nWARNING: isotope Z=" << z() << ", A=" << a()
              << " - cannot find rrms. Using approx formula: rrms="
              << approx_rrms << "\n\n";
  }

  if (t() < 0.0) {
    t() = default_t;
  }
}

Nucleus::Nucleus(const std::string &z_str, int in_a,
                 const std::string &str_type, double in_rrms, double in_t,
                 const std::vector<double> &in_params,
                 const std::string & custom_pot_file_name)
    : Nucleus(AtomData::atomic_Z(z_str), in_a, str_type, in_rrms, in_t,
              in_params, custom_pot_file_name) {}

//==============================================================================
std::ostream &operator<<(std::ostream &ostr, const Nucleus &n) {
  const auto rrms = n.r_rms();
  const auto tt = n.t();

  switch (n.type()) {
  case Nuclear::ChargeDistro::point:
    ostr << "Point-like nucleus; ";
    break;
  case Nuclear::ChargeDistro::spherical:
    ostr << "Spherical nucleus; "
         << " r_rms = " << rrms
         << ", r_charge = " << Nuclear::c_hdr_formula_rrms_t(rrms, 0);
    break;
  case Nuclear::ChargeDistro::Gaussian:
    ostr << "Gaussian nucleus; "
         << " r_rms = " << rrms
         << ", r_charge = " << Nuclear::c_hdr_formula_rrms_t(rrms, 0);
    break;
  case Nuclear::ChargeDistro::Fermi: {
    const auto c = Nuclear::c_hdr_formula_rrms_t(rrms, tt);
    ostr << "Fermi nucleus; "
         << " r_rms = " << rrms
         << ", c_hdr = " << c
         << ", t = " << tt;

        if(n.beta() != 0.0) {
          ostr << ", \u03b2 = " << n.beta()
               << ", t_eff = "  << Nuclear::deformation_effective_t(c, tt, n.beta());
        }
    break;
  }
    case Nuclear::ChargeDistro::custom:
    ostr << "Custom nuclear potential; "
         << "loaded from file: \'" << n.custom_pot_file()
         << "\'";
    break;
  default:
    ostr << "Invalid nucleus type?";
  }

  if (!n.params().empty()) {
    ostr << " (parameters: ";
    for (auto p : n.params())
      ostr << p << ", ";
    ostr << ")";
  }

  return ostr;
}

//==============================================================================
Nucleus form_nucleus(int Z, std::optional<int> tA, IO::InputBlock input) {

  // Nucleus: Get + setup nuclear parameters

  input.check(
      {{"",   "Options for nuclear potential (finite nuclear size). All "
            "are optional. Default is a Fermi-like nucleus, with "
            "parameters chosen according to isotope (see Atom{A;})"},
       {"rrms", "Root-mean-square charge radius, in fm. If not given, will use "
            "default for given Z and A."},
       {"c", "Half-density radius, in fm. ( If c and rrms are defined, c will "
            "over-ride rms). Default depends on Z and A."},
       {"type", "Fermi, spherical, point-like, Gaussian [Fermi], custom"},
       {"t", "Nuclear skin thickness, in fm [2.3]"},
       {"beta", "Intrinsic nuclear quadrupole deformation parameter \u03b2 "
            "(Currently only for Fermi potential) [0.0]"},
       {"input_file", "Input text file containing nuclear potential grid values to use."},
       {"parameters",
        "List of comma separated real numbers. Not currently unused, but may "
        "be used in future for more complicated nuclear distros."}});

  using namespace std::string_literals;
  // Set nuclear type if given
  const auto nuc_type = input.get("type", "Fermi"s);

  // Use A. If not given, find default for Z
  const auto A = tA ? *tA : AtomData::defaultA(Z);

  // Get default nucleus:
  auto nucleus = Nuclear::Nucleus{Z, A, nuc_type};

  // over-ride default options if given
  const auto rrms = input.get<double>("rrms");
  const auto t = input.get<double>("t");
  const auto c_hdr = input.get<double>("c");
  const auto beta = input.get<double>("beta");

  const auto input_file = input.get<std::string>("input_file");

  if (t) {
    nucleus.t() = *t;
  }
  if (rrms) {
    nucleus.set_rrms(*rrms);
  }
  if (c_hdr) {
    // this will over-ride given rms
    nucleus.set_rrms(Nuclear::rrms_formula_c_t(*c_hdr, nucleus.t()));
  }
  if (beta) {
      nucleus.beta() = *beta;
  }

  if (input_file) {
    nucleus.custom_pot_file() = *input_file;
  }

  // Set "extra" parameters (not currently used)
  nucleus.params() = input.get<std::vector<double>>("parameters", {});

  // If A or given rrms are zero, explicitly set to point-like nucleus
  // This isn't required, but makes output more explicit
  if (nucleus.a() == 0.0 || nucleus.r_rms() == 0.0) {
    nucleus.t() = 0.0;
    nucleus.set_rrms(0.0);
    nucleus.type() = Nuclear::ChargeDistro::point;
  }
  return nucleus;
}

//==============================================================================
std::vector<double> sphericalNuclearPotential(double Z, double rnuc,
                                              const std::vector<double> &rgrid)
// Potential due to a spherical nucleus, with (charge) radius, rnuc.
// Note: rnuc must be given in "fermi" (fm, femto-metres).
// rnuc = 0 corresponds to zeroNucleus
{
  std::vector<double> vnuc;
  vnuc.reserve(rgrid.size());

  // Fill the vnuc array with spherical nuclear potential
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

//==============================================================================
std::vector<double> GaussianNuclearPotential(double Z, double r_rms,
                                             const std::vector<double> &rgrid)

{
  std::vector<double> vnuc;
  vnuc.reserve(rgrid.size());

  // Fill the vnuc array with Gaussian nuclear potential
  const double rN = r_rms / PhysConst::aB_fm; // convert fm -> au
  const auto k = std::sqrt(3.0 / 2);
  for (auto r : rgrid) {
    vnuc.push_back(-Z * std::erf(k * r / rN) / r);
  }

  return vnuc;
}

//==============================================================================
std::vector<double> fermiNuclearPotential(double z, double t, double c,
                                          const std::vector<double> &rgrid)
// Uses a Fermi-Dirac distribution for the nuclear potential.
//
// rho(r) = rho_0 {1 + Exf[(r-c)/a]}^-1
// V(r) = -(4 Pi)/r [A+B]
//   A = Int[ rho(x) x^2 , {x,0,r}]
//   B = r * Int[ rho(x) x , {x,r,infty}]
// rho_0 is found by either:
//   * V(infinity) = -Z/r , or equivalently
//   * int rho(r) d^3r = Z
//
// Depends on:
//   * t: skin thickness [90 to 10% fall-off range]
//     note: t = a[4 ln(3)]
//   * c: half-density radius (see above for link c and r_rms)
//
// t and c are input values. In 'fermi' or fm (femto metres)
//
// V(r) is expressed in terms of Complete Fermi-Dirac integrals.
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
    const double roc = PhysConst::aB_fm * r / c;
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

//==============================================================================
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

  using namespace qip::overloads;
  rho *= rho0;

  return rho;
}

std::vector<double> readCustomNuclearPotential(const std::string& if_name, std::size_t rgrid_size)
{
  std::vector<double> vnuc;
  vnuc.reserve(rgrid_size);
  
  std::ifstream potential_fstream;
  potential_fstream.open(if_name);
  double current_value = 0.0;

  // Consistency check to ensure the potential grid has the same number of points as the radial grid.
  std::size_t number_lines = 0;
  while (potential_fstream >> current_value) {
      vnuc.push_back(current_value);
      ++number_lines;
  }

  if(number_lines != rgrid_size) {
    std::cout << "Error: provided nuclear potential file: \"" << if_name << "\" has incompatible size"
              << " (expected " << rgrid_size << ", got " << number_lines << ").";
    std::abort();
  }

  return vnuc;
}

//==============================================================================
std::vector<double> formPotential(const Nucleus &nuc,
                                  const std::vector<double> &r) {
  const auto z = nuc.z();
  const auto r_rms = nuc.r_rms();
  const auto t = nuc.t();

  // const auto params = nuc.params(); // not used as yet

  const auto nucleus_type =
      r_rms > 0.0 ? nuc.type() : Nuclear::ChargeDistro::point;

  switch (nucleus_type) {

  case Nuclear::ChargeDistro::Fermi: {
    const auto chdr = Nuclear::c_hdr_formula_rrms_t(r_rms, t);
    const auto t_eff = (nuc.beta() != 0.0) ? deformation_effective_t(chdr, t, nuc.beta()) : t;
    return Nuclear::fermiNuclearPotential(z, t_eff, chdr, r);
  }

  case Nuclear::ChargeDistro::spherical: {
    const auto r_n = std::sqrt(5.0 / 3.0) * r_rms;
    return Nuclear::sphericalNuclearPotential(z, r_n, r);
  }

  case Nuclear::ChargeDistro::point: {
    // return Nuclear::sphericalNuclearPotential(z, 0.0, r);
    using namespace qip::overloads;
    return -z / r;
  }

  case Nuclear::ChargeDistro::Gaussian: {
    return Nuclear::GaussianNuclearPotential(z, r_rms, r);
  }

  case Nuclear::ChargeDistro::custom: {
    return Nuclear::readCustomNuclearPotential(nuc.custom_pot_file(), r.size());
  }

  default: {
    std::cout << "\nFAIL Nuclear:266 - invalid nucleus type?\n";
    std::abort();
  }
  }
}

} // namespace Nuclear