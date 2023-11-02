#pragma once
#include "NuclearData.hpp" //for Isotope
#include <string>
#include <vector>
class Grid;

namespace Nuclear {

//! Nuclear charge distribution options
enum class ChargeDistro { Fermi, spherical, point, Gaussian };

ChargeDistro parseType(const std::string &str_type);
std::string parseType(ChargeDistro type);

//------------------------------------------------------------------------------

//! Stores set of nuclear parameters (all radii in fm)
class Nucleus {
  Isotope m_iso;
  Nuclear::ChargeDistro m_type;
  double m_t;

public:
  Nucleus(int in_z = 1, int in_a = 0, const std::string &str_type = "Fermi",
          double in_rrms = -1.0, double in_t = -1.0);

  Nucleus(const std::string &z_str, int in_a,
          const std::string &str_type = "Fermi", double in_rrms = -1.0,
          double in_t = Nuclear::default_t);

public:
  ChargeDistro &type() { return m_type; }
  ChargeDistro type() const { return m_type; }
  int &z() { return m_iso.Z; };
  int z() const { return m_iso.Z; };
  int &a() { return m_iso.A; };
  int a() const { return m_iso.A; };
  // std::optional<double> &r_rms() { return m_iso.r_rms; };

  void set_rrms(double rrms) { m_iso.r_rms = rrms; }

  double r_rms() const { return m_iso.r_rms ? *m_iso.r_rms : 0.0; };
  double &t() { return m_t; };
  double t() const { return m_t; };

  double c() const { return c_hdr_formula_rrms_t(r_rms(), m_t); }

  friend std::ostream &operator<<(std::ostream &ostr, const Nucleus &n);
};

//------------------------------------------------------------------------------

//! Nuclear potentials: spherical charge distribution
[[nodiscard]] std::vector<double>
sphericalNuclearPotential(double Z, double rnuc,
                          const std::vector<double> &rgrid);

//! Nuclear potentials: Gaussian charge distribution
[[nodiscard]] std::vector<double>
GaussianNuclearPotential(double Z, double r_rms,
                         const std::vector<double> &rgrid);

//! Nuclear potentials: Fermi charge distribution [c is half-density radius, not
//! rms]
[[nodiscard]] std::vector<double>
fermiNuclearPotential(double Z, double t, double c,
                      const std::vector<double> &rgrid);

//! Fermi charge distribution, rho(r) - normalised to Z_norm
[[nodiscard]] std::vector<double>
fermiNuclearDensity_tcN(double t, double c, double Z_norm, const Grid &grid);

//! Calls one of the above, depending on params. Fills V(r), given r
[[nodiscard]] std::vector<double> formPotential(const Nucleus &nucleus,
                                                const std::vector<double> &r);

} // namespace Nuclear
