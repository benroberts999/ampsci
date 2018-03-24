#include "PRM_parametricPotentials.h"


//******************************************************************************
double PRM_green(int Z, double r, double H, double d)
/*
Have default H and d : scaled to fit Alkali atoms!
Does not include nuclear part!
*/
{
  if(r==0) return 0;
  return (-1./r)*(1. + (double(Z)-1.)/(H*(exp(r/d)-1)+1)) + double(Z)/r;
}


//******************************************************************************
double PRM_tietz(int Z, double r, double t, double g)
/*
TIETZ parametric potential
Does NOT include nuclear potential
*/
{
  if(r==0) return 0;
	return (-1./r)*(1. + (double(Z)-1.)*exp(-g*r)/pow(1.+t*r,2)) + double(Z)/r;
}



//******************************************************************************
int PRM_defaultGreen(int z, double &H, double &d)
/*
Default values for Green potential parameters.
Determined by fitting the 5 lowest J=1/2 states.
Crude quadratic fit used for other Z values.
*/
{

  if(z==3){
    H=0.65074;
    d=0.37017;
    return 0;
  }else if(z==11){
    H=1.32493;
    d=0.48946;
    return 0;
  }else if(z==19){
    H=2.49872;
    d=0.78867;
    return 0;
  }else if(z==37){
    H=3.60424;
    d=0.80155;
    //Johnson:
    // H=3.4811;
    // d=0.7855;
    return 0;
  }else if(z==49){
    H=3.16518;
    d=0.71266;
    return 0;
  }else if(z==55){
    H=4.83108;
    d=0.93920;
    //Johnson:
    // H=4.4691;
    // d=0.8967;
    return 0;
  }else if(z==79){
    //Johnson:
    H=4.4560;
    d=0.7160;
    return 0;
  }else if(z==81){
    H=4.44226;
    d=0.72244;
    //Johnson:
    // H=4.4530;
    // d=0.7234;
    return 0;
  }else if(z==87){
    H=7.02125;
    d=1.00633;
    return 0;
  }

  //If not one of the defaults, use crude model fit
  H = 0.665356 + 0.0747058*z - 0.000155671*pow(z,2);
  d = 0.433183 + 0.0103617*z - 0.0000531958*pow(z,2);
  return 0;
}

//******************************************************************************
int PRM_defaultTietz(int z, double &t, double &g)
/*
Default values for Green potential parameters.
Determined by fitting the 5 lowest J=1/2 states.
Crude quadratic fit used for other Z values.
*/
{
  if(z==3){
    t=0.10000;
    g=1.94363;
    return 0;
  }else if(z==5){
    t=0.81084;
    g=0.24086;
    return 0;
  }else if(z==11){
    t=0.60635;
    g=1.38417;
    return 0;
  }else if(z==13){
    t=1.22936;
    g=0.26552;
    return 0;
  }else if(z==19){
    t=1.13831;
    g=0.55698;
    return 0;
  }else if(z==31){
    t=2.05009;
    g=0.18498;
    return 0;
  }else if(z==37){
    t=1.65346;
    g=0.46553;
    // Johnson:
    // t=1.9530;
    // g=0.2700;
    return 0;
  }else if(z==49){
    t=2.36372;
    g=0.09911;
    return 0;
  }else if(z==55){
    t=1.91561;
    g=0.31948;
    // Johnson:
    // t=2.0453;
    // g=0.2445;
    return 0;
  }else if(z==79){
    t=2.4310;
    g=0.3500;
    return 0;
  }else if(z==81){
    // t=2.90583;
    // g=0.10749;
    // Johnson:
    t=2.3537;
    g=0.3895;
    return 0;
  }else if(z==87){
    t=2.03525;
    g=0.50664;
    return 0;
  }

  //If not one of the defaults, use crude model fit
  t = 0.33991 + 0.0495393*z - 0.000272772*pow(z,2);
  g = 1.14316 - 0.028863*z + 0.000205713*pow(z,2);
  return 0;
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
