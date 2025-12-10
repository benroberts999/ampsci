#pragma once
#include "Physics/PhysConst_constants.hpp"

namespace UnitConv {

//! energy: au -> keV
constexpr double Energy_au_to_keV = PhysConst::Hartree_eV / 1.0e3;
//! energy: keV -> au
constexpr double Energy_keV_to_au = 1.0 / Energy_au_to_keV;

//! energy: au -> keV
constexpr double Energy_au_to_eV = PhysConst::Hartree_eV;
//! energy: keV -> au
constexpr double Energy_eV_to_au = 1.0 / Energy_au_to_eV;

//! energy: cm^-1 -> au
constexpr double Energy_invcm_to_au = 1.0 / PhysConst::Hartree_invcm;

//! momentum: au -> MeV: [hbar*q] = [hbar/a0] = (m_e*c*alpha) = E_H/c*alpha
constexpr double Momentum_au_to_MeV =
    PhysConst::Hartree_eV / PhysConst::alpha / 1.0e6;
//! momentum: MeV -> au
constexpr double Momentum_MeV_to_au = 1.0 / Momentum_au_to_MeV;

//! momentum: au -> eV: [hbar*q] = [hbar/a0] = (m_e*c*alpha) = E_H/c*alpha
constexpr double Momentum_au_to_eV = PhysConst::Hartree_eV / PhysConst::alpha;
//! momentum: eV -> au
constexpr double Momentum_eV_to_au = 1.0 / Momentum_au_to_eV;

//! mass: au -> GeV
constexpr double Mass_au_to_GeV = PhysConst::m_e_MeV / 1000.0;
//! mass: au -> MeV
constexpr double Mass_au_to_MeV = PhysConst::m_e_MeV;
//! mass: GeV -> au
constexpr double Mass_GeV_to_au = 1.0 / Mass_au_to_GeV;
//! mass: MeV -> au
constexpr double Mass_MeV_to_au = 1.0 / Mass_au_to_MeV;

//! velocity: au -> km/s
constexpr double Velocity_au_to_kms = (PhysConst::c_SI / PhysConst::c) / 1.0e3;
//! velocity: au -> cm/s
constexpr double Velocity_au_to_cms =
    (PhysConst::c_SI / PhysConst::alpha) * 100.;
//! velocity: cm/s -> cm/day
constexpr double Velocity_au_to_cmday =
    Velocity_au_to_cms * (24.0 * 60.0 * 60.0);

//! atomic mass: daltons -> kg
constexpr double AtomicMass_u_to_kg = PhysConst::u_NMU * PhysConst::m_e_kg;

} // namespace UnitConv