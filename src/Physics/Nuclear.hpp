#pragma once
#include "Nuclear_DataTable.hpp" //for Isotope
#include <string>
#include <vector>
class Grid;

namespace Nuclear {

// skin-thickness. Always same?
constexpr double default_t = 2.30;

enum class Type { Fermi, spherical, point };

Type parseType(const std::string &str_type);

struct Parameters {
  int z, a;
  Nuclear::Type type;
  double t;
  double r_rms;
  Parameters(const std::string &z_str, int in_a,
             const std::string &str_type = "Fermi", double in_rrms = -1.0,
             double in_t = -1.0);
  Parameters(int in_z, int in_a, const std::string &str_type = "Fermi",
             double in_rrms = -1.0, double in_t = -1.0);
};

//******************************************************************************
Isotope findIsotopeData(int z, int a);
std::vector<Isotope> findIsotopeList(int z);
double find_rrms(int z, int a);
double find_mu(int z, int a);
int find_parity(int z, int a);
double find_spin(int z, int a);

//******************************************************************************
double approximate_r_rms(int a);
double c_hdr_formula_rrms_t(double rrms, double t = default_t);
double rrms_formula_c_t(double c, double t = default_t);
constexpr double approximate_t_skin(int a);

std::vector<double> sphericalNuclearPotential(double Z, double rnuc,
                                              const std::vector<double> &rgrid);

std::vector<double> fermiNuclearPotential(double Z, double t, double c,
                                          const std::vector<double> &rgrid);

std::vector<double> fermiNuclearDensity_tcN(double t, double c, double Z_norm,
                                            const Grid &grid);

std::vector<double> formPotential(Parameters params, int z, int,
                                  const std::vector<double> &r);

} // namespace Nuclear
