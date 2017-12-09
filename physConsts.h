#ifndef _PHYSCONSTS_H
#define _PHYSCONSTS_H
//#include <cmath>


double const CLIGHT=137.035999074;				//speed of light in a.u.
double const CLIGHT2=CLIGHT*CLIGHT;
double const ALPHA=1./CLIGHT;						// fine structure constant
double const ALPHA2=ALPHA*ALPHA;
//// kill below!
double const cc=137.035999074;				//speed of light in a.u.
double const c2=cc*cc;
double const aa=1./cc;						// fine structure constant
double const aa2=aa*aa;

// units conversions:
double const ABOHR=0.52917721092e-10;			// Bohr radius (in metres) 	(1 au = ab m)
double const RYDBERG=13.60569253;				// 1 Rydberg (in eV) (1H=2Ry)

double const aB=0.52917721092e-10;			// Bohr radius (in metres) 	(1 au = ab m)
double const Ryd=13.60569253;				// 1 Rydberg (in eV) (1H=2Ry)
double const efund=1.602176565e-19;			// fundamental charge C, [1eV=efund J]
double const ccsi=299792458.;				// s.o.l. in SI m/s
double const hbarc=197.3269718;				// hbar*c, in MeV.fm
double const invcm=2.*Ryd*8065.54429;		// Inverse cm:	(1 au = invcm cm^-1) [1H * eV(in 1/cm)]

// electron mass
double const meMeV=0.510998928;				// mass_electron in MeV/c^2
double const mesi=9.10938291e-31;			// " "  in SI (kg)

#endif
