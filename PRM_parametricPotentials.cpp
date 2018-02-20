#include "PRM_parametricPotentials.h"


//******************************************************************************
double PRM_green(int Z, double r, double H, double d)
/*
Have default H and d : scaled to fit Alkali atoms!
Does not include nuclear part!
*/
{
	//double iH=4.4691;		//params for Cs
	//double id=0.8967;
//	double iH=3.4811;		//params for Rb
//	double id=0.7855;
  if(r==0) return 0;
  return (-1./r)*(1. + (double(Z)-1.)/(H*(exp(r/d)-1)+1)) + double(Z)/r;

}


//******************************************************************************
double PRM_tietz(int Z, double r, double g, double t)
/*
TIETZ parametric potential
Does NOT include nuclear potential
*/
{
	// double ig=0.2445;		//params for Cs
	// double it=2.0453;
  if(r==0) return 0;
	return (-1./r)*(1. + (double(Z)-1.)*exp(-g*r)/pow(1.+t*r,2)) + double(Z)/r;
}

// //******************************************************************************
// // "Approximate" (expansion) solution for Thomas-Fermi potential for neutral atom
// double TFapprox(int i, int Z){
// // for pure coulomb.. - infinitesimal nucleus (can change)
// // Agrees very will with GREEN parametric potential [Checked for Cs only]
// // Potential corresponds to V^N-1 potential
// 	double x;
// 	double phi;
// 	double vTF;
// 	if (i==0){
// //		vTF=-Z/(r0);			// ?? or let it be infinite??.. This shouldn't be used
// 		vTF=-Z/(0.01*r0);			// ?? or let it be infinite??.. This shouldn't be used
// 	}
// 	else if ((i>0)and(i<NGP)){
// //		x=1.12950781018323*r(i)*pow(Z,1/3.);
// 		x=1.12950781018323*r(i)*cbrt(Z);
// 		phi=1./(1+0.02747*pow(x,0.5)+1.243*x-0.1486*pow(x,1.5)
//           +0.2302*pow(x,2)+0.007298*pow(x,2.5)+0.006944*pow(x,3)
//         );
// //		vTF = -((Z-1)*phi+1)/r(i);
// 		vTF = -(0.95)*((Z-1)*phi+1)/r(i);
// 	}
// 	else{
// 		printf("FAILURE: Wrong gridpoint used for TFapprox(i=%i)\n",i);
// 		return 0;
// 	}
// 	return vTF;
// }
