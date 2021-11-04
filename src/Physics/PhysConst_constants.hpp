#pragma once
#include <cmath> // PI

//! @brief Set of commonly-used Physics constants
//! @details
//! Constains physical constants, and units conversions.
//! Taken mostly from 2018 CODATA values.
//! https://physics.nist.gov/cuu/Constants/
namespace PhysConst {

// CODATA 2018: 1/alpha = 137.035999084(21)
//! speed of light in a.u.
constexpr double c = 137.035999084;
constexpr double c2 = c * c;
//! fine-structure constant
constexpr double alpha = 1.0 / c;
constexpr double alpha2 = alpha * alpha;

// speed of light, in m/s [CODATA 2018, exact]
constexpr double c_SI = 299792458.0;
constexpr double hbar_SI = (6.62607015e-34) / (2 * M_PI); // exact

//! CODATA 2018: hbar*c = 197.326 980 4... MeV fm
constexpr double hbarc_MeVfm = 197.3269804;

// Mass (mp/me)
// CODATA 2018: m_e = 9.109 383 7015(28) e-31
// CODATA 2018: m_p = 1.672 621 923 69(51) x 10-27 kg
// codata 2018: 1836.152 673 43(11)
//! Proton mass, in atomic units (mp/me)
constexpr double m_p = 1836.15267343;
constexpr double m_e_kg = 9.1093837015e-31;
//! MeV/c^2 - electron mass; 2020 value 0.51099895000(15):
constexpr double m_e_MeV = 0.51099895000;

// CODATA 2014: 1822.888 486 192(53)
//! "unified atomic mass" unit; nuclear mass unit; Dalton; u
constexpr double u_NMU = 1822.888486192;

// Length:
// CODATA 2018: a_B = 0.529177210903(80)e-10 m
constexpr double aB_m = 0.529177210903e-10;
constexpr double aB_cm = aB_m * (1.0e+2);
constexpr double aB_fm = aB_m * (1.0e+15);

// Time: hbar/E_H = wiki: 2.4188843265857(47)×10−17 s [wiki 2019]
//! Time: hbar/E_H
constexpr double time_s = 2.4188843265857e-17;

// Energy:
// CODATA 2018: E_H = 27.211 386 245 988(53) eV
//                  = 6.579 683 920 502(13) x 10^15 Hz
constexpr double Hartree_eV = 27.211386245988;
constexpr double Hartree_Hz = 6.579683920502e+15;
constexpr double Hartree_MHz = Hartree_Hz * (1.0e-6);
constexpr double Hartree_GHz = Hartree_Hz * (1.0e-9);
// wave-number (inverse cm):
constexpr double Hartree_invcm = 1.0 / (2.0 * M_PI * c * aB_cm);
// wavelength (nm):
constexpr double HartreeWL_nm = 2.0 * M_PI * c * aB_m * (1.0e+9);

//! Fermi weak constant (au). Times 10^{11}
/*! Particle Data Group 2020:
 Gf = 1.1663787(6) e-5 (hbar*c)^3 GeV^-2
    = 1.1663787(6) e-5 alpha (GeV/c^2)^-2 au
 me = 0.51099895000 e-3 GeV/c^2 = 1 au
 Gf = 2.222516(11) e-3 au
*/
constexpr double GFe11 = 2.222516e-3;
constexpr double GF = GFe11 * (1.0e-11);

//! electron g-factor. CODATA 2018: -2.002 319 304 362 56(35)
constexpr double g_e = -2.00231930436256;

//! Bohr magneton (in SI-derived atomic units):
constexpr double muB_SI = 0.5; //
//! Bohr magneton (in Gaussian CGS-derived atomic units):
constexpr double muB_CGS = 0.5 * alpha;
//! Nulcear magneton (in SI-derived atomic units):
constexpr double muN_SI = muB_SI / m_p;
//! Nulcear magneton (in Gaussian CGS-derived atomic units):
constexpr double muN_CGS = muB_CGS / m_p;
//! Nulcear magneton in MHz (via Gaussian CGS-derived atomic units):
constexpr double muN_CGS_MHz = Hartree_MHz * muB_CGS / m_p;

//! barn = 1.0e-28m^2, for Quadrupole moment
constexpr double barn_m2 = 1.0e-28;
constexpr double barn_au = barn_m2 / (aB_m * aB_m);
constexpr double barn_MHz = barn_au * Hartree_MHz;

} // namespace PhysConst
