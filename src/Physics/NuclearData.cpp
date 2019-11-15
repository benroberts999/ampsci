#include "Physics/NuclearData.hpp" //for Isotope

namespace Nuclear {
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
  return (A < 10) ? 1.15 * std::pow(A, 0.3333)
                  : 0.836 * std::pow(A, 0.3333) + 0.570;
}

//******************************************************************************
double c_hdr_formula_rrms_t(double rrms, double t)
// Calculates half-density radius, given rms charge radius, and t.
// r_rms = sqrt(<r^2>) = (3/5)c^2 + (7/5)(pi*a)^2
// a = t / (4 * ln[3])
{
  const double a = t / FourLn3;
  return (rrms < t)
             ? std::sqrt((5.0 / 3.0) * rrms * rrms)
             // this is little dodgy? but formula prob  only works large A
             : std::sqrt((5.0 / 3.0) * rrms * rrms -
                         (7.0 / 3.0) * (Pi2 * a * a));
}

//******************************************************************************
double rrms_formula_c_t(double c, double t)
// Calculates  rms charge radius, given half-density radius (c), and t.
// r_rms = sqrt(<r^2>) = (3/5)c^2 + (7/5)(pi*a)^2
// a = t / (4 * ln[3])
{
  const double a = t / FourLn3;
  return std::sqrt((3.0 / 5.0) * c * c + (7.0 / 5.0) * a * a * Pi2);
}

//******************************************************************************

double approximate_t_skin(int) { return default_t; }

} // namespace Nuclear
