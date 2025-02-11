#pragma once
#include "IO/InputBlock.hpp"
#include "Physics/NuclearData.hpp" //for Isotope
#include <optional>
#include <string>
#include <vector>
class Grid;

namespace Nuclear {

//! Nuclear charge distribution options
enum class ChargeDistro { Fermi, spherical, point, Gaussian, custom, Error };

ChargeDistro parseType(const std::string &str_type);
std::string parseType(ChargeDistro type);

//------------------------------------------------------------------------------

//! Stores set of nuclear parameters (all radii in fm)
class Nucleus {
  // Isotope data
  Isotope m_iso;
  // Charge distribution type
  Nuclear::ChargeDistro m_type;
  // skin thickness parameter (in fm). Usually 2.3
  double m_t;
  // Deformation parameter beta. Usually 0.
  double m_beta;
  // Name of input file used to read in custom nuclear potential
  std::string m_custom_pot_file_name;
  // Other parameters: not used for now
  std::vector<double> m_params;

public:
  Nucleus(int in_z = 1, int in_a = 0, const std::string &str_type = "Fermi",
          double in_rrms = -1.0, double in_t = -1.0, double in_beta = 0.0,
          const std::vector<double> &in_params = {},
          const std::string &custom_pot_file_name = "");

  Nucleus(const std::string &z_str, int in_a,
          const std::string &str_type = "Fermi", double in_rrms = -1.0,
          double in_t = Nuclear::default_t, double in_beta = 0.0,
          const std::vector<double> &in_params = {},
          const std::string &custom_pot_file_name = "");

public:
  ChargeDistro &type() { return m_type; }
  ChargeDistro type() const { return m_type; }
  int &z() { return m_iso.Z; };
  int z() const { return m_iso.Z; };
  int &a() { return m_iso.A; };
  int a() const { return m_iso.A; };

  const std::vector<double> &params() const { return m_params; };
  std::vector<double> &params() { return m_params; };

  void set_rrms(double rrms) { m_iso.r_rms = rrms; }
  double r_rms() const { return m_iso.r_rms ? *m_iso.r_rms : 0.0; };

  double &t() { return m_t; };
  double t() const { return m_t; };

  double &beta() { return m_beta; };
  double beta() const { return m_beta; };

  std::string &custom_pot_file() { return m_custom_pot_file_name; };
  std::string custom_pot_file() const { return m_custom_pot_file_name; };

  double c() const { return c_hdr_formula_rrms_t(r_rms(), m_t); }

  friend std::ostream &operator<<(std::ostream &ostr, const Nucleus &n);
};

//------------------------------------------------------------------------------
Nucleus form_nucleus(int Z, std::optional<int> A = std::nullopt,
                     IO::InputBlock input = {});

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
