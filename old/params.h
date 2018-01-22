#ifndef _PARAMS_H
#define _PARAMS_H
#include <cmath>
/******************************************
Wed 05 Nov 2014 13:09:25 AEDT

Global Parameters for all programs.

*******************************************/
/*
double r0=0.0005; 			// (0.0005)
int NGP=500;      		//probably THIS should be a function of r_max and h and r0 etc.  (500)
double h=0.03125;			//this is the one that should be as small as possible! (0.03125)
*/

//params for radial grid!
//const double paramRmax=50000;  // this should be INPUT?
//const double r0=5e-6;
//const double h=0.005;			// for some reason, this doesn't work if h<0.00048.. ?
//const int NGP=log(paramRmax/r0)/h+2;

//const double paramRmax=100;  // this should be INPUT?
const double paramRmax=500;  // this should be INPUT?
const double r0=1e-4;
const int NGP=2000;
const double h=log(paramRmax/r0)/(NGP-2);

// NOTE maxing these out, doesn't increase accuracy that much.... and breaks other things, such as finite nuclear size!
// e.g. 0.01 gives same accuracy as 0.0005 !
// but (for up to n=25 (ALL kappa) for H, with r0=1e-6:
// 0.01 takes  0.701658 s  [r0=1e-3 :: 0.483439]
// 0.0005 takes 13.305596 seconds  [r0=1e-3 :: 8.938231 (one didn't conv.)]

// bound state wavefunctions (Adams-moul)
const double delep=5e-15;		//PRIMARY convergence parameter for bound state energy	(10^-11)
const double deles=1e-11;		//SECONDAY convergence parameter for bound state energy	(X)
const int ntry=30;			// Number of failed attempts at converging (sove bs) before error quit (30)
const double alr=800;			// ''assymptotically large r [not what this is..]''  (=800)
const double lde=0.2;		//amount to vary energy by for 'large' variations (0.1 => 10%)
// (in-dir)
const int nx=30;				//order of the expansion coeficients in 'inint'  (15 orig.)
const float nxepsp=1e-18;		// PRIMARY convergance for expansion in `inint'' (10^-8)
const float nxepss=1e-3;		// SECONDARY convergance for expansion in `inint'' (10e-3)
// (out-dir)
const int nol=1;				// # of outdir runs [finds first nol*AMO+1 points (3)]

const int nquad=14; 	// any number from 1-14. for 'nquad'-point quadrature integration


//int AMO=8; // this should be a const.... but then forces 'seg fault' onto ie[][].. away from coefa.....???????
const int amo2=7;		// When this is 8, it doesn't work well !!! ??
const int AMO=amo2;

/////////////////////////////////////////////////////////////////
// DATA:

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
