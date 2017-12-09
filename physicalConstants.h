#ifndef _PHYSCONSTS_H
#define _PHYSCONSTS_H
/*
Constains physical constants, and units conversions.
Taken mostly from 2014 CODATA values.
https://physics.nist.gov/cuu/Constants/
https://en.wikipedia.org/wiki/Atomic_units
*/

//speed of light in a.u., and fine-structure constant
double const CLIGHT=137.035999139;  //CODATA 2014: 137.035999139(31)
double const CLIGHT2=CLIGHT*CLIGHT;
double const ALPHA=1./CLIGHT; // fine structure constant
double const ALPHA2=ALPHA*ALPHA;

//Proton mass (mp/me)
const double MPROTON=1836.15267389; //CODATA 2014: 1836.152 673 89(17)

//Length:
double const ABOHR_M=0.52917721067e-10; //CODATA 2014: 0.52917721067(12)e-10 m
double const ABOHR_CM=0.52917721067e-8;

//Time:
double const TIME_S=2.418884326505e-17; //wiki: 2.418884326505(16)×10−17 s

//Energy:
double const HARTREE_EV=27.21138602; //CODATA 2014: 27.21138602(17) eV
double const HARTREE_HZ=6.579683920711e+15; // 6.579683920711(39)e15 Hz
double const HARTREE_MHZ=6.579683920711e+9;
double const HARTREE_GHZ=6.579683920711e+6;
// wave-number (inverse cm):
double const HARTREE_ICM=2.194746313702e+5; //2.194746313702(13)e7 m-1
// wavelength (nm):
double const HARTREE_WLNM=45.56335252767;

// Bohr magneton (in atomic units):
const double MUB_SI=0.5;        //SI-derived
const double MUB_CGS=0.5*ALPHA; //Gaussian CGS-derived
// Nulcear magneton (in atomic units):
const double MUN_SI=MUB_SI/MPROTON;
const double MUN_CGS=MUB_CGS/MPROTON;

#endif
