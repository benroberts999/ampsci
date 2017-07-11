/******************************************
Wed 05 Nov 2014 13:09:25 AEDT 

Global `Data' for all programs.
-Includes fundamental constants, and untis conversions

********************************************/
// !! Make cc able to be varied!

double const cc=137.035999074;				//speed of light in a.u.
double const c2=cc*cc;
double const aa=1/cc;						// fine structure constant
double const aa2=aa*aa;

// units conversions:
double const aB=0.52917721092e-10;			// Bohr radius (in metres) 	(1 au = ab m)
double const Ryd=13.60569253;				// 1 Rydberg (in eV) (1H=2Ry)
double const efund=1.602176565e-19;			// fundamental charge C, [1eV=efund J]
double const ccsi=299792458;				// s.o.l. in SI m/s
double const hbarc=197.3269718;				// hbar*c, in MeV.fm
double const invcm=2*Ryd*8065.54429;		// Inverse cm:	(1 au = invcm cm^-1) [1H * eV(in 1/cm)]

// electron mass
double const meMeV=0.510998928;				// mass_electron in MeV/c^2
double const mesi=9.10938291e-31;			// " "  in SI (kg)
