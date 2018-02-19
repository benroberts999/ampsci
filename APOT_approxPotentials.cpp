#include "APOT_approxPotentials.h"


//******************************************************************************
double APOT_green(int Z, double r, double H, double d)
/*
Have default H and d : scaled to fit Alkali atoms!
XXX do not include nuclear part!
*/
{
// for pure coulomb.. - infinitesimal nucleus
	double vgreen;
	//double iH=4.4691;		//params for Cs
	//double id=0.8967;
//	double iH=3.4811;		//params for Rb
//	double id=0.7855;
  if(r==0) return 0;
  return (-1./r)*(1. + (double(Z)-1.)/(H*(exp(r/d)-1)+1)) + double(Z)/r;

}
