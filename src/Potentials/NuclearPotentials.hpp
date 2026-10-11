#pragma once
#include "IO/InputBlock.hpp"
#include "Physics/NuclearData.hpp" //for Isotope
#include <optional>
#include <string>
#include <vector>
class Grid;

namespace Nuclear {

/*!
  @brief Nuclear charge distribution options.

  @details
  - Fermi     : two-parameter Fermi distribution (c and t)
  - spherical : uniformly charged sphere (ball)
  - point     : point-like nucleus
  - Gaussian  : Gaussian distribution
  - custom    : potential read from file (see formPotential)
  - Error     : unrecognised type (see parseType)
*/
enum class ChargeDistro { Fermi, spherical, point, Gaussian, custom, Error };

/*!
  @brief Converts string to ChargeDistro (case-insensitive).

  @details
  Also accepts: "spher*", "ball", "uniform" (spherical); "point*", "zero*"
  (point); "gaus*" (Gaussian). If unknown, prints a warning with the closest
  match.

  @return ChargeDistro; ChargeDistro::Error if unknown.
*/
ChargeDistro parseType(const std::string &str_type);
//! Converts ChargeDistro to string
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
  /*!
    @brief Constructs nucleus from Z and A, with optional parameters.

    @details
    r_rms is taken from isotope data if not given; if not in the isotope
    data, an approximate formula is used. A = 0 gives r_rms = 0 (point-like),
    regardless of @p in_rrms. For non-Fermi types, t is set to zero.

    @param in_z                  Nuclear charge, Z
    @param in_a                  Mass number, A. Negative: default A for Z;
                                 0: point-like
    @param str_type              Charge distribution type (see parseType)
    @param in_rrms               Root-mean-square charge radius (fm).
                                 Negative: from isotope data
    @param in_t                  Skin thickness (fm), Fermi only. Negative:
                                 default_t
    @param in_beta               Quadrupole deformation parameter (Fermi only)
    @param in_params             Extra parameters (not currently used)
    @param custom_pot_file_name  File for custom nuclear potential
  */
  Nucleus(int in_z = 1, int in_a = 0, const std::string &str_type = "Fermi",
          double in_rrms = -1.0, double in_t = -1.0, double in_beta = 0.0,
          const std::vector<double> &in_params = {},
          const std::string &custom_pot_file_name = "");

  //! As above, with Z given as atomic symbol or name (e.g., "Cs")
  Nucleus(const std::string &z_str, int in_a,
          const std::string &str_type = "Fermi", double in_rrms = -1.0,
          double in_t = Nuclear::default_t, double in_beta = 0.0,
          const std::vector<double> &in_params = {},
          const std::string &custom_pot_file_name = "");

public:
  //! Charge distribution type
  ChargeDistro &type() { return m_type; }
  ChargeDistro type() const { return m_type; }
  //! Nuclear charge, Z
  int &z() { return m_iso.Z; };
  int z() const { return m_iso.Z; };
  //! Mass number, A
  int &a() { return m_iso.A; };
  int a() const { return m_iso.A; };

  //! Extra parameters (not currently used)
  const std::vector<double> &params() const { return m_params; };
  std::vector<double> &params() { return m_params; };

  //! Sets root-mean-square charge radius (fm)
  void set_rrms(double rrms) { m_iso.r_rms = rrms; }
  //! Root-mean-square charge radius (fm); 0 if not set
  double r_rms() const { return m_iso.r_rms ? *m_iso.r_rms : 0.0; };

  //! Skin thickness, t (fm); Fermi only
  double &t() { return m_t; };
  double t() const { return m_t; };

  //! Quadrupole deformation parameter, beta; Fermi only
  double &beta() { return m_beta; };
  double beta() const { return m_beta; };

  //! File name for custom nuclear potential
  std::string &custom_pot_file() { return m_custom_pot_file_name; };
  std::string custom_pot_file() const { return m_custom_pot_file_name; };

  //! Half-density radius, c (fm), from r_rms and t
  double c() const { return c_hdr_formula_rrms_t(r_rms(), m_t); }

  //! Prints summary of nuclear parameters
  friend std::ostream &operator<<(std::ostream &ostr, const Nucleus &n);
};

//------------------------------------------------------------------------------
/*!
  @brief Forms nucleus from Z, A, and the nuclear options in @p input.

  @details
  Starts from the default nucleus for given Z and A (see Nucleus), then
  over-rides with any options given in @p input: type, rrms, c, t, beta,
  input_file, parameters. If c is given, it over-rides rrms. If A or r_rms
  is zero, the nucleus is set to point-like (except for custom).

  @param Z      Nuclear charge
  @param A      Mass number; default A for Z if not given
  @param input  Input block with nuclear options (see input.check() in .cpp)

  @note input_file only sets the file name: type must also be custom.
*/
Nucleus form_nucleus(int Z, std::optional<int> A = std::nullopt,
                     IO::InputBlock input = {});

//------------------------------------------------------------------------------

/*!
  @brief Nuclear potential: uniformly charged sphere (ball) of radius rnuc.

  @details
  \f[
    V(r) = \begin{cases}
      -\frac{Z}{2 r_N}\left(3 - r^2/r_N^2\right) & r < r_N, \\
      -Z/r & r \geq r_N.
    \end{cases}
  \f]
  rnuc = 0 gives a point-like nucleus.

  @param Z      Nuclear charge
  @param rnuc   Nuclear charge radius, in fm (not r_rms: for a uniform
                sphere, r_N = sqrt(5/3) r_rms)
  @param rgrid  Radial points (atomic units)

  @return V(r) at each point of @p rgrid (atomic units).
*/
[[nodiscard]] std::vector<double>
sphericalNuclearPotential(double Z, double rnuc,
                          const std::vector<double> &rgrid);

/*!
  @brief Nuclear potential: Gaussian charge distribution.

  @details
  \f[
    V(r) = -\frac{Z}{r}\,{\rm erf}\left(\sqrt{3/2}\;r/r_{\rm rms}\right).
  \f]

  @param Z      Nuclear charge
  @param r_rms  Root-mean-square charge radius, in fm
  @param rgrid  Radial points (atomic units)

  @return V(r) at each point of @p rgrid (atomic units).
*/
[[nodiscard]] std::vector<double>
GaussianNuclearPotential(double Z, double r_rms,
                         const std::vector<double> &rgrid);

/*!
  @brief Nuclear potential: Fermi charge distribution.

  @details
  For charge density
  \f[
    \rho(r) = \rho_0 \left[ 1 + e^{(r - c)/a} \right]^{-1},
    \qquad a = t/(4\ln 3),
  \f]
  normalised to Z. Expressed in terms of complete Fermi-Dirac integrals
  (GSL). V = -Z/r beyond r = 10c.

  @param Z      Nuclear charge
  @param t      Skin thickness (90% to 10% fall-off), in fm
  @param c      Half-density radius, in fm (not r_rms)
  @param rgrid  Radial points (atomic units)

  @return V(r) at each point of @p rgrid (atomic units).
*/
[[nodiscard]] std::vector<double>
fermiNuclearPotential(double Z, double t, double c,
                      const std::vector<double> &rgrid);

/*!
  @brief Fermi charge distribution, rho(r), normalised to Z_norm.

  @details
  \f$ \int \rho\, d^3r = 4\pi \int \rho(r)\, r^2 dr = Z_{\rm norm} \f$:
  Z_norm = Z for nuclear charge density; Z_norm = 1 for nuclear density.
  rho = 0 beyond r = 10c.

  @param t       Skin thickness, in fm
  @param c       Half-density radius, in fm
  @param Z_norm  Normalisation
  @param grid    Radial grid

  @return rho(r) on @p grid (atomic units).
*/
[[nodiscard]] std::vector<double>
fermiNuclearDensity_tcN(double t, double c, double Z_norm, const Grid &grid);

/*!
  @brief Nuclear potential V(r) for given nucleus, at radial points @p r.

  @details
  Calls one of the above, depending on nucleus type. If r_rms = 0, point-like
  (-Z/r) for all types except custom.
  - Fermi: c from r_rms and t; if beta is non-zero, uses the effective skin
    thickness for the deformed nucleus (deformation_effective_t)
  - spherical: radius r_N = sqrt(5/3) r_rms
  - custom: reads {r, V(r)} (atomic units) from custom_pot_file(), and
    interpolates onto @p r. Beyond the last point, V = -Z_eff/r, with Z_eff
    from the last point. Not extrapolated to smaller r.

  @return V(r) at each point of @p r (atomic units).

  @warning Aborts if the custom file is invalid (fewer than 3 points), or the
           type is invalid.
*/
[[nodiscard]] std::vector<double> formPotential(const Nucleus &nucleus,
                                                const std::vector<double> &r);

} // namespace Nuclear
