#ifndef _PHYSCONSTS_H
#define _PHYSCONSTS_H
/*
Constains physical constants, and units conversions.
Taken mostly from 2014 CODATA values.
https://physics.nist.gov/cuu/Constants/
https://en.wikipedia.org/wiki/Atomic_units
*/

namespace FPC{

//speed of light in a.u., and fine-structure constant
double const c=137.035999139;  //CODATA 2014: 137.035999139(31)
double const c2=c*c;
double const alpha=1./c; // fine structure constant
double const alpha2=alpha*alpha;

//Proton mass (mp/me)
const double m_p=1836.15267389; //CODATA 2014: 1836.152 673 89(17)
const double m_e_MeV = 0.5109989461; //MeV/c^2 - electron mass

//Length:
double const aB_m=0.52917721067e-10; //CODATA 2014: 0.52917721067(12)e-10 m
double const aB_cm=0.52917721067e-8;
double const aB_fm=0.52917721067e+5;

//Time:
double const time_s=2.418884326505e-17; //wiki: 2.418884326505(16)×10−17 s

//Energy:
double const Hartree_eV=27.21138602; //CODATA 2014: 27.21138602(17) eV
double const Hartree_Hz=6.579683920711e+15; // 6.579683920711(39)e15 Hz
double const Hartree_MHz=6.579683920711e+9;
double const Hartree_GHz=6.579683920711e+6;
// wave-number (inverse cm):
double const Hartree_invcm=2.194746313702e+5; //2.194746313702(13)e7 m-1
// wavelength (nm):
double const HartreeWL_nm=45.56335252767;

// Bohr magneton (in atomic units):
const double muB_SI=0.5;        //SI-derived
const double muB_CGS=0.5*alpha; //Gaussian CGS-derived
// Nulcear magneton (in atomic units):
const double muN_SI=muB_SI/m_p;
const double muN_CGS=muB_CGS/m_p;

}
#endif
