#pragma once
#include <cmath>
#include <optional>
#include <vector>

//! Data and useful functions for nuclear properties and potentials. Radii all
//! in Fermi (fm, e-15m) from Nuclear Data Service: https://www-nds.iaea.org/
namespace Nuclear {

// skin-thickness. Always same?
constexpr double default_t = 2.30;

constexpr auto FourLn3 = 4.0 * 1.098612289;
// 4.0 * std::log(3.0); // std::log not constepr
constexpr auto Pi2 = M_PI * M_PI;

//! Isotope data: Z, A, r_rms/fm, I, pi, mu, Q
struct Isotope {
  //! Atomic charge
  int Z;
  //! Atomic mass number (A = Z + N)
  int A;
  //! root-mean-square charge radius, in Fermi (fm, e-15m)
  std::optional<double> r_rms{};
  //! Nuclear spin (in hbar)
  std::optional<double> I_N{};
  std::optional<int> parity{};
  //! Magnetic dipole moment, in nuclear magnetons
  std::optional<double> mu{};
  //! Magnetic quadrupole moment, in barns
  std::optional<double> q{};
};

//==============================================================================
//! Looks up + returns an isotope from the list. If not in list, partially blank
Isotope findIsotopeData(int z, int a);
//! Returns all known isotopes of given atom
std::vector<Isotope> findIsotopeList(int z);
//! Looks up default value of r_rms for given isotope. Returns 0 if not found.
double find_rrms(int z, int a);
//! As above, for dipole moment
double find_mu(int z, int a);
//! As above, for parity
int find_parity(int z, int a);
//! As above, for nuclear spin. Returns -1 if not found
double find_spin(int z, int a);

//==============================================================================
//! Calculates c from rrms and t
double c_hdr_formula_rrms_t(double rrms, double t = default_t);
//! Calculates rrms from c and t
double rrms_formula_c_t(double c, double t = default_t);

//! Aproximate rms radius from a fir to Angeli data
double approximate_r_rms(int a, int z);
//! just returns 2.3
double approximate_t_skin(int a);

} // namespace Nuclear
