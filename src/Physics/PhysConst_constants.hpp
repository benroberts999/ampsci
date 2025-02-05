#pragma once
#include <cmath> // PI

//! @brief Set of commonly-used Physics constants
//! @details
//! Constains physical constants, and units conversions.
//! Taken mostly from 2018 CODATA values.
//! https://physics.nist.gov/cuu/Constants/
namespace PhysConst {

//! c, speed of light: 299 792 458 m/s   [exact]
constexpr double c_SI = 299792458.0;

//! Planck constant h: 6.626 070 15 e-34 J.s [exact]
constexpr double h_SI = (6.62607015e-34);

//! hbar: 6.626 070 15 e-34 / (2 Pi) J.s [exact]
constexpr double hbar_SI = h_SI / (2.0 * M_PI);

//! Fundamental charge, Coulombs. 1.602 176 634 e-19 C [exact]
constexpr double e_C = 1.602176634e-19;

//! Fine-structure constant: alpha = 1/137.035 999 177(21) [CODATA 2022]
constexpr double alpha = 1.0 / 137.035999177;

//! electron g-factor: -2.002 319 304 360 92(36) [CODATA 2022]
constexpr double g_e = -2.00231930436092;

//! e/me: 1.758 820 008 38(55) e11 [CODATA 2022]
constexpr double e_on_me_SI = 1.75882000838e11;

//! electron mass, in SI (kg)
constexpr double m_e_kg = e_C / e_on_me_SI;

//! Proton mass, in atomic units (mp/me). CODATA 2022: 1836.152 673 426(32)
constexpr double m_p = 1836.152673426;

//! muon mass in atomic units (m_mu/m_e). Codata 2022: 206.768 2827(46)
constexpr double m_muon = 206.7682827;

//! tauon mass in atomic units (m_tau/m_e) 3477.23(23)
constexpr double m_tau = 3477.23;

//! unified atomic mass unit; (nuclear mass unit, Dalton): CODATA 2022: 1.660 539 068 92(52) x 10-27 kg
constexpr double u_NMU_kg = 1.66053906892e-27;

//! Bohr radius, in m. : hbar/(m_e*c*alpha)
constexpr double aB_m = hbar_SI / (m_e_kg * c_SI * alpha);

//! Hartree (atomic energy unit = 2Ry) in eV, 27.211 386 245 981(30) eV [CODATA 2022]
constexpr double Hartree_eV = 27.211386245981;
// nb: this one different from m*c^2*alpha^2 - probably because it measured better than alpha??

//! Fermi weak constant, in GeV^-2: 1.166 3787(6) x 10-5 GeV-2 [CODATA 2022]
constexpr double GF_GeV2 = 1.1663787e-5;

//======================================================================
//======================================================================

constexpr double alpha2 = alpha * alpha;

//! speed of light in a.u. (=1/alpha)
constexpr double c = 1.0 / alpha;
constexpr double c2 = c * c;

//! hbar * c, in MeV.fm
constexpr double hbarc_MeVfm = (hbar_SI * c_SI / e_C) * 1.0e9;

//! Electron mass (MeV/c^2)
constexpr double m_e_MeV = Hartree_eV * c2 / 1.0e6;

//! unified atomic mass unit; (nuclear mass unit, Dalton): au
constexpr double u_NMU = u_NMU_kg / m_e_kg;

// Length:

constexpr double aB_cm = aB_m * (1.0e+2);
constexpr double aB_fm = aB_m * (1.0e+15);
constexpr double aB_nm = aB_m * (1.0e+9);

//! Hartree (atomic energy unit = 2Ry) in Hz.
constexpr double Hartree_Hz = Hartree_eV * e_C / h_SI;

constexpr double Hartree_MHz = Hartree_Hz * (1.0e-6);
constexpr double Hartree_GHz = Hartree_Hz * (1.0e-9);
//! Hartree to cm^-1 conversion [wave-number, inverse cm]:
constexpr double Hartree_invcm = 1.0 / (2.0 * M_PI * c * aB_cm);
//! Hartree to corresponding wavelength, in nm
constexpr double HartreeWL_nm = 2.0 * M_PI * c * aB_m * (1.0e+9);

//! hbar/E_H (atomic unit of time) (in seconds)
constexpr double hbar_on_EH = hbar_SI / e_C / Hartree_eV;

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

//! Fermi weak constant (au).
/*! Particle Data Group 2020:
 Gf = 1.1663787(6) e-5 (hbar*c)^3 GeV^-2
    = 1.1663787(6) e-5 alpha (GeV/c^2)^-2 au
 me = 0.51099895000 e-3 GeV/c^2 = 1 au
 Gf = 2.222516(11) e-3 au
*/
constexpr double GF = GF_GeV2 * alpha * m_e_MeV * m_e_MeV * 1e-6;
//! Fermi weak constant * 10^11, in atomic units
constexpr double GFe11 = GF * 1.0e11;

} // namespace PhysConst
