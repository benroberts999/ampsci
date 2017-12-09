#ifndef _PHYSCONSTS_H
#define _PHYSCONSTS_H
//#include <cmath>


//speed of light in a.u.
double const CLIGHT=137.035999139;  //CODATA 2014: 137.035999139(31)
double const CLIGHT2=CLIGHT*CLIGHT;
double const ALPHA=1./CLIGHT;						// fine structure constant
double const ALPHA2=ALPHA*ALPHA;

// units conversions (from atomic to _):

//Length:
double const ABOHR_SI=0.52917721067e-10; //CODATA 2014: 0.52917721067(12)e-10 m

//Energy:
double const HARTREE_EV=27.21138602; //CODATA 2014: 27.21138602(17) eV
double const HARTREE_HZ=6.579683920711e+15; // 6.579683920711(39)e15 Hz
// wave-number (inverse cm):
double const HARTREE_ICM=2.194746313702e+5; //2.194746313702(13)e7 m-1
// wavelength (nm):
double const HARTREE_WLNM=45.56335252767;


#endif
