#include "INT_quadratureIntegration.h"
#include <vector>
//*********************************************************
//*********************************************************
//*********************************************************

namespace INT{


//Coefs for 14-point int
int nquad14 = 14;
double c14[14] =
{1360737653653,9821965479386,-7321658717812,34616887868158,
-48598072507095,81634716670404,-78837462715392,76782233435964,-41474518178601,
28198302087170,-3009613761932,7268021504806,4920175305323,5252701747968};
double dd14 = 5230697472000;



//			INTEGRATE
// Integrates input function f(i) wrt r, from l to m. This program contains wronskian drdt
//Uses an nquad-point quadrature formula.. for nquad=1->14 (any integer)
//double integrate(double *f, int l, int m)
double integrate(const std::vector<double> &fin, const std::vector<double> &w, double h, int l,
  int m, int nquad)
/*

f(x) dx -> f(u(t)) w(t) dt

w = du/dt = dr/dt is the wronskian

XXX overload so w is optional!
XXX overload so can use floats?

*/
{

//Note: normally, NGP is so large, that there is no difference between nquad=1 to 14!

  std::vector<double>f=fin; //XXX bad temp. workaround! XXX

  int ngp = f.size();
  if (m==0) m=ngp-1;

//re-arange limits to make less fiddly
// so I don't need to ..
	int negInt=1;
	if (l>m){
		int temp=l;
		l=m;
		m=temp;
		negInt=-1;
	}
	if (m>=ngp){m=ngp-1;}
	if (l<0){l=0;}


  if(nquad<1)nquad=1;
  if(nquad>14)nquad=14;

// Defines the `nquad'-point quadrature integration coeficents.
// nquad can take any integer from 1 to 14 (=1 implies trapazoid rule)
	//double ic[nquad];
	std::vector<double> ic(nquad);
	if (nquad==1){
		double c[1]={1};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==2){
		double c[2]={1,2};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==3){
		double c[3]={9,28,23};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==4){
		double c[4]={8,31,20,25};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==5){
		double c[5]={475,1902,1104,1586,1413};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==6){
		double c[6]={459,1982,944,1746,1333,1456};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==7){
		double c[7]={36799,176648,54851,177984,89437,130936,119585};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==8){
		double c[8]={35584,185153,29336,220509,46912,156451,111080,122175};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==9){
		double c[9]={2082753,11532470,261166,16263486,-1020160,12489922,5095890,
		7783754,7200319};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==10){
		double c[10]={2034625,11965622,-1471442,20306238,-7084288,18554050,1053138,
		9516362,6767167,7305728};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==11){
		double c[11]={262747265,1637546484,-454944189,3373884696,-2145575886,3897945600,
		-1065220914,1942518504,636547389,1021256716,952327935};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==12){
		double c[12]={257696640,1693103359,-732728564,4207237821,-3812282136,6231334350,
		-3398609664,3609224754,-196805736,1299041091,896771060,963053825};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==13){
		double c[13]={1382741929621,9535909891802,-5605325192308,28323664941310,
		-32865015189975,53315213499588,-41078125154304,39022895874876,-13155015007785,
		12465244770050,3283609164916,5551687979302,5206230892907};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==14){
		double c[14]={1360737653653,9821965479386,-7321658717812,34616887868158,
		-48598072507095,81634716670404,-78837462715392,76782233435964,-41474518178601,
		28198302087170,-3009613761932,7268021504806,4920175305323,5252701747968};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else {
		//printf("FAILURE: Wrong order for integration.. check nquad\n");
		return 1; //XXX
	}

	double dd[14]={2,2,24,24,1440,1440,120960,120960,7257600,
					7257600,958003200,958003200,5230697472000,5230697472000};


	// wronskian
  //XXX make optional!
	for (int i=l; i<=m; i++){
		f[i]=f[i]*w[i]; //XXX double check - is this OK??
	}

	// Defines/calculates the integral
	double a=0;
	int ll=l;
	int mm=m;
	for (int i=0; i<nquad; i++){
		a=a+ic[i]*(f[ll]+f[mm]);
		ll=ll+1;
		mm=mm-1;
	}
	a=a/dd[nquad-1];
	for (int i=l; i<=m; i++){
		a=a+f[i];
	}
	a=a*h*negInt;

	// // rectangle method (for tests)
	// double b=0;
	// for (int i=l; i<m; i++){
	// 	b=b+f[i];
	// }
	// b=b*h;

	return a;
}		// END integrate





//******************************************************************
// Calculates the derivative of a function!
int diff(const std::vector<double> &f, const std::vector<double> &drdt, double h,
    std::vector<double> &deriv)
{

// Problem: This calcs the derivative for a few (4) points AFTER pinf; where it oscilates.
// Means the first 0 of dP won't correspond to first max of P!

// This hasn't been fully tested!

/*
df/dr = df/dt * dt/dr = (df/dt) / (dr/dt) = (df/di) * (di/dt) / (dr/dt) =
= (df/di)  / (h * dr/dt)
coeficients from: http://en.wikipedia.org/wiki/Finite_difference_coefficient
*/

  int NGP = f.size();

	deriv[0]= (f[1]-f[0])/(h*drdt[0]);
	deriv[NGP-1]= (f[NGP-1]-f[NGP-2])/(h*drdt[NGP-1]);

	deriv[1]= (f[2]-f[0])/(2*h*drdt[1]);
	deriv[NGP-2]= (f[NGP-1]-f[NGP-3])/(2*h*drdt[NGP-2]);

	deriv[2]= (f[0]-8*f[1]+8*f[3]-f[4])/(12*h*drdt[2]);
	deriv[NGP-3]= (f[NGP-5]-8*f[NGP-4]+8*f[NGP-2]-f[NGP-1])/(12*h*drdt[NGP-3]);

	deriv[3]= (-1*f[0]+9*f[1]-45*f[2]+45*f[4]-9*f[5]+1*f[6])/(60*h*drdt[3]);
	deriv[NGP-4]= (-1*f[NGP-7]+9*f[NGP-6]-45*f[NGP-5]+45*f[NGP-3]-9*f[NGP-2]+1*f[NGP-1])/(60*h*drdt[NGP-4]);

	for (int i=4; i<(NGP-4);i++){
//		deriv[i] = ( (1/56)*f[i-4] - (4/21)*f[i-3] + 1*f[i-2] - 4*f[i-1]
//		- (1/56)*f[i+4] + (4/21)*f[i+3] - 1*f[i+2] + 4*f[i+1] ) / (5*h*drdt(i));
		deriv[i] = ( (1./8)*f[i-4] - (4./3)*f[i-3] + 7*f[i-2] - 28*f[i-1]
		- (1./8)*f[i+4] + (4./3)*f[i+3] - 7*f[i+2] + 28*f[i+1] ) / (35*h*drdt[i]);
	}

  return 0;

}

//******************************************************************
double integrate1(
  const std::vector<double> &f1,
  double h,
  int l, int m,
  int a_start, int a_end)
/*
Integrates f1*f2*f3
NOTE: f3 should be wronskian!Jacobian?Whatever
*/
{
  int ngp = f1.size();
  if(m==0) m = ngp-1;
  if (m>=ngp)  m = ngp-1;
  if (m<=l) return 0;
	// Defines/calculates the integral
	double a = 0;
  for (int i=l; i<=m; i++) a += f1[i];

  if( (a_start==0 && a_end==0) || ((m-l) <= 2*nquad14)) return a*h;

	int ll = l;
	int mm = m ;
	for (int i=0; i<nquad14; i++){
		a += (c14[i]/dd14)*(a_start*f1[ll] + a_end*f1[mm]);
		ll++;
		mm--;
	}

	return a*h;
}		// END integrate1
//******************************************************************
double integrate2(
  const std::vector<double> &f1,
  const std::vector<double> &f2,
  double h,
  int l, int m,
  int a_start, int a_end)
/*
Integrates f1*f2*f3
NOTE: f3 should be wronskian!Jacobian?Whatever
*/
{
  int ngp = f1.size();
  if(m==0) m = ngp-1;
  if (m>=ngp)  m = ngp-1;
  if (m<=l) return 0;
	// Defines/calculates the integral
  double a = 0;
  for (int i=l; i<=m; i++) a += f1[i]*f2[i];

  if( (a_start==0 && a_end==0) || ((m-l) <= 2*nquad14)) return a*h;

	int ll = l;
	int mm = m - 1;
	for (int i=0; i<nquad14; i++){
		a += (c14[i]/dd14)*(a_start*f1[ll]*f2[ll] + a_end*f1[mm]*f2[mm]);
		ll++;
		mm--;
	}

	return a*h;
}		// END integrate2
//******************************************************************
double integrate3(
  const std::vector<double> &f1,
  const std::vector<double> &f2,
  const std::vector<double> &f3,
  double h,
  int l, int m,
  int a_start, int a_end)
/*
Integrates f1*f2*f3
NOTE: f3 should be wronskian!Jacobian?Whatever
*/
{
  int ngp = f1.size();
  if(m==0) m = ngp-1;
  if (m>=ngp)  m = ngp-1;
  if (m<=l) return 0;
	// Defines/calculates the integral
  double a = 0;
  for (int i=l; i<=m; i++) a += f1[i]*f2[i]*f3[i];

  if( (a_start==0 && a_end==0) || ((m-l) <= 2*nquad14)) return a*h;

	int ll = l;
	int mm = m ;
	for (int i=0; i<nquad14; i++){
		a += (c14[i]/dd14)*(a_start*f1[ll]*f2[ll]*f3[ll]
        + a_end*f1[mm]*f2[mm]*f3[mm]);
		ll++;
		mm--;
	}

	return a*h;
}		// END integrate3
//******************************************************************
double integrate4(
  const std::vector<double> &f1,
  const std::vector<double> &f2,
  const std::vector<double> &f3,
  const std::vector<double> &f4,
  double h,
  int l, int m,
  int a_start, int a_end)
/*
Integrates f1*f2*f3
NOTE: f3 should be wronskian!Jacobian?Whatever
*/
{
  int ngp = f1.size();
  if(m==0) m = ngp-1;
  if (m>=ngp)  m = ngp-1;
  if (m<=l) return 0;
  // Defines/calculates the integral
  double a = 0;
  for (int i=l; i<=m; i++) a += f1[i]*f2[i]*f3[i]*f4[i];

  if( (a_start==0 && a_end==0) || ((m-l) <= 2*nquad14)) return a*h;

	int ll = l;
	int mm = m ;
	for (int i=0; i<nquad14; i++){
		a += (c14[i]/dd14)*(a_start*f1[ll]*f2[ll]*f3[ll]*f4[ll]
        + a_end*f1[mm]*f2[mm]*f3[mm]*f4[mm]);
		ll++;
		mm--;
	}

	return a*h;
}		// END integrate3







}//namespace
