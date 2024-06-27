#include "Physics/NuclearData.hpp" //for Isotope
#include "Physics/nuclear_data_table.hpp"

namespace Nuclear {
//==============================================================================
Isotope findIsotopeData(int z, int a) {
  if (a < z)
    return Isotope{z, 0, 0.0, 0, {}, {}, {}};
  for (const auto &nucleus : NuclearDataTable) {
    if (nucleus.Z == z && nucleus.A == a)
      return nucleus;
  }
  return Isotope{z, a, {}, {}, {}, {}, {}};
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
  const auto nuc = findIsotopeData(z, a);
  return nuc.r_rms ? *nuc.r_rms : 0.0;
}

double find_mu(int z, int a) {
  auto nuc = findIsotopeData(z, a);
  return nuc.mu ? *nuc.mu : 0.0;
}

int find_parity(int z, int a) {
  const auto nuc = findIsotopeData(z, a);
  return nuc.parity ? *nuc.parity : 0;
}

double find_spin(int z, int a) {
  const auto nuc = findIsotopeData(z, a);
  return nuc.I_N ? *nuc.I_N : 0.0;
}

//==============================================================================
double approximate_r_rms(int A, int Z)
// Returns approximate root-mean-square charge radius in fm [1.e-15 m]
// https://www.sciencedirect.com/science/article/pii/S0092640X12000265
// https://www-nds.iaea.org/radii/
{
  // return (A < 10) ? 1.15 * std::pow(A, 0.3333) :
  //                   0.836 * std::pow(A, 0.3333) + 0.570;

  const auto a13 = std::pow(A, 1.0 / 3.0);
  const auto z13 = std::pow(Z, 1.0 / 3.0);

  // From mathematica fit to Angeli data
  if (A <= 0)
    return 0.0;
  double rt{0.0};
  if (A % 2 == 0) {
    rt = (A < 20) ? 1.15875 - 0.384013 * a13 + 1.26460 * z13 :
                    0.344936 + 0.554395 * a13 + 0.439873 * z13;
  } else {
    rt = (A < 5) ? -1.90973 + 1.99164 * a13 + 0.796396 * z13 :
                   0.488201 + 0.60582 * a13 + 0.333936 * z13;
  }
  // round to 2 digits so as not to give impression of accuracy
  return (static_cast<int>(rt * 100.0)) / 100.0;
}

//==============================================================================
double c_hdr_formula_rrms_t(double rrms, double t)
// Calculates half-density radius, given rms charge radius, and t.
// r_rms = sqrt(<r^2>) = (3/5)c^2 + (7/5)(pi*a)^2
// a = t / (4 * ln[3])
{
  const double a = t / FourLn3;
  return (rrms < t) ?
             std::sqrt(5.0 / 3.0) * rrms
             // this is little dodgy? but formula prob  only works large A
             :
             std::sqrt((5.0 / 3.0) * rrms * rrms - (7.0 / 3.0) * (Pi2 * a * a));
}

//==============================================================================
double rrms_formula_c_t(double c, double t)
// Calculates  rms charge radius, given half-density radius (c), and t.
// r_rms = sqrt(<r^2>) = (3/5)c^2 + (7/5)(pi*a)^2
// a = t / (4 * ln[3])
{
  const double a = t / FourLn3;
  return std::sqrt((3.0 / 5.0) * c * c + (7.0 / 5.0) * a * a * Pi2);
}

double deformation_effective_t(double c, double t, double beta)
// Nuclear quadrupole deformation is approximately equivalent
// to change in skin thickness (see, e.g., Eq. 8 10.1103/PhysRevA.100.032511)
{
  const double sqrt_delta_t = FourLn3 * c * beta;
  return std::sqrt(t * t + 3.0 / (4.0 * M_PI * M_PI * M_PI) * sqrt_delta_t * sqrt_delta_t);
}

//==============================================================================

double approximate_t_skin(int) { return default_t; }

} // namespace Nuclear
