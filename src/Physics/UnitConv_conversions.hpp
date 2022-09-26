#pragma once
#include "Physics/PhysConst_constants.hpp"

namespace UnitConv {

// energy: au -> keV
constexpr double Energy_au_to_keV = PhysConst::Hartree_eV / 1.0e3;
// energy: keV -> au
constexpr double Energy_keV_to_au = 1. / Energy_au_to_keV;
//! energy: cm^-1 -> au
constexpr double Energy_invcm_to_au = 1. / PhysConst::Hartree_invcm;

//! momentum: au -> MeV
constexpr double Momentum_au_to_MeV =
    PhysConst::Hartree_eV * PhysConst::c / 1.0e6;
//! momentum: MeV -> au
constexpr double Momentum_MeV_to_au = 1. / Momentum_au_to_MeV;

//! mass: au -> GeV
constexpr double Mass_au_to_GeV = PhysConst::m_e_MeV / 1000.;
//! mass: au -> MeV
constexpr double Mass_au_to_MeV = PhysConst::m_e_MeV;
//! mass: GeV -> au
constexpr double Mass_GeV_to_au = 1. / Mass_au_to_GeV;
//! mass: MeV -> au
constexpr double Mass_MeV_to_au = 1. / Mass_au_to_MeV;

//! velocity: au -> km/s
constexpr double Velocity_au_to_kms = (PhysConst::c_SI / PhysConst::c) / 1.0e3;
//! velocity: au -> cm/s
constexpr double Velocity_au_to_cms =
    (PhysConst::c_SI / PhysConst::alpha) * 100.;
//! velocity: cm/s -> cm/day
constexpr double Velocity_au_to_cmday = Velocity_au_to_cms * (24. * 60. * 60.);

//! atomic mass: daltons -> kg
constexpr double AtomicMass_u_to_kg = PhysConst::u_NMU * PhysConst::m_e_kg;

} // namespace UnitConv