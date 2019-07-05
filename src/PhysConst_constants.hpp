#pragma once
#include <cmath> // PI
/*
Constains physical constants, and units conversions.
Taken mostly from 2018 CODATA values.
https://physics.nist.gov/cuu/Constants/
https://en.wikipedia.org/wiki/Atomic_units
*/

namespace PhysConst {

// speed of light in a.u., and fine-structure constant
// CODATA 2018: 1/alpha = 137.035999084(21)
const double c = 137.035999084;
const double c2 = c * c;
const double alpha = 1. / c; // fine structure constant
const double alpha2 = alpha * alpha;

// speed of light, in m/s [CODATA 2018, exact]
const double c_SI = 299792458.;
const double hbar_SI = (6.62607015e-34) / (2 * M_PI); // exact

// Mass (mp/me)
// CODATA 2018: m_e = 9.109 383 7015(28) e-31
// CODATA 2018: m_p = 1.672 621 923 69(51) x 10-27 kg
const double m_p = 1836.15267343; // codata 2018: 1836.152 673 43(11)
const double m_e_kg = 9.1093837015e-31;
const double m_e_MeV = 0.5109989461; // MeV/c^2 - electron mass; 2014 value

//"unified atomic mass" unit; nuclear mass unit; Dalton; u
const double u_NMU = 1822.888486192; // CODATA 2014: 1822.888 486 192(53)

// Length:
// CODATA 2018: a_B = 0.529177210903(80)e-10 m
const double aB_m = 0.529177210903e-10;
const double aB_cm = aB_m * (1.0e+2);
const double aB_fm = aB_m * (1.0e+15);

// Time: hbar/E_H = wiki: 2.4188843265857(47)×10−17 s [wiki 2019]
const double time_s = 2.4188843265857e-17;

// Energy:
// CODATA 2018: E_H = 27.211 386 245 988(53) eV
//                  = 6.579 683 920 502(13) x 10^15 Hz
const double Hartree_eV = 27.211386245988;
const double Hartree_Hz = 6.579683920502e+15;
const double Hartree_MHz = Hartree_Hz * (1.0e-6);
const double Hartree_GHz = Hartree_Hz * (1.0e-9);
// wave-number (inverse cm):
const double Hartree_invcm = 1.0 / (2.0 * M_PI * c * aB_m * (1.0e+2));
// wavelength (nm):
const double HartreeWL_nm = 2.0 * M_PI * c * aB_m * (1.0e+9);

// Fermi weak constant (au)
const double GFe11 = 2.2225e-3;
const double GF = GFe11 * (1.0e-11);

// electron g-factor
// CODATA 2018: -2.002 319 304 362 56(35)
const double g_e = -2.00231930436256;

// Bohr magneton (in atomic units):
const double muB_SI = 0.5;          // SI-derived
const double muB_CGS = 0.5 * alpha; // Gaussian CGS-derived
// Nulcear magneton (in atomic units):
const double muN_SI = muB_SI / m_p;
const double muN_CGS = muB_CGS / m_p;

} // namespace PhysConst
