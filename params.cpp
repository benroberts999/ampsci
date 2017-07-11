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
const double paramRmax=50000;  // this should be INPUT?
const double r0=5e-6;
const double h=0.005;			// for some reason, this doesn't work if h<0.00048.. ?
const int NGP=log(paramRmax/r0)/h+2;

// NOTE maxing these out, doesn't increase accuracy that much.... and breaks other things, such as finite nuclear size!
// e.g. 0.01 gives same accuracy as 0.0005 !
// but (for up to n=25 (ALL kappa) for H, with r0=1e-6:
// 0.01 takes  0.701658 s  [r0=1e-3 :: 0.483439]
// 0.0005 takes 13.305596 seconds  [r0=1e-3 :: 8.938231 (one didn't conv.)]



// bound state wavefunctions (Adams-moul)
double delep=5e-15;		//PRIMARY convergence parameter for bound state energy	(10^-11)
double deles=1e-11;		//SECONDAY convergence parameter for bound state energy	(X)
int ntry=30;			// Number of failed attempts at converging (sove bs) before error quit (30)
double alr=800;			// ''assymptotically large r [not what this is..]''  (=800)
double lde=0.2;		//amount to vary energy by for 'large' variations (0.1 => 10%)
// (in-dir)
int nx=30;				//order of the expansion coeficients in 'inint'  (15 orig.)
float nxepsp=1e-18;		// PRIMARY convergance for expansion in `inint'' (10^-8)
float nxepss=1e-3;		// SECONDARY convergance for expansion in `inint'' (10e-3)
// (out-dir)
int nol=1;				// # of outdir runs [finds first nol*AMO+1 points (3)]

const int nquad=14; 	// any number from 1-14. for 'nquad'-point quadrature integration


//int AMO=8; // this should be a const.... but then forces 'seg fault' onto ie[][].. away from coefa.....???????
const int amo2=7;		// When this is 8, it doesn't work well !!! ??
int AMO=amo2;








