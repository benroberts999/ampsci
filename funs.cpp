//////////////////////////////////////////////////////////////////////

/*
********************************************************
Tue 20 Jan 2015  

Contains functions that can be called from other programs:

-radial gridpoint

-nuclear and parametric potentials

-integration

- Exact (0-size nucleus) Dirac energies, for comparison

-An "approximation" to the Thomas-Fermi potential (based on an expansion of Z-dependent solution. 


*********************************************************
*/



//*********************************************************
//*********************************************************
//*********************************************************
//			DEFINES RADIAL GRID

// dr/dt
double drdt(int i){			
	double rpi;
	if (i==0){
		rpi=r0;
	}
	else if ((i>0)and(i<NGP)){
		rpi=r0*exp((i)*h);
	}
	else{
		printf("FAILURE: Wrong gridpoint used for drdt(i=%i)\n",i);
		return 0;
	}
	return rpi;
}

// defines radial grid
double r(int i){
	double ri;
	if (i==0){
		ri=0.001*r0;
	}
	else if ((i>0)and(i<NGP)){
		ri=drdt(i)-r0;
	}
	else{
		printf("FAILURE: Wrong gridpoint used for r(i=%i)\n",i);
		return 0;
	}
	return ri;
}

// (dr/dt)/r
double dror(int i){
	double rpori;
	if (i==0){
		rpori=0;
	}
	else if ((i>0)and(i<NGP)){
		rpori=drdt(i)/r(i);
	}
	else{
		printf("FAILURE: Wrong gridpoint used for dror(i=%i)\n",i);
		return 0;
	}
	return rpori;
}








//*********************************************************
//*********************************************************
//*********************************************************
//			POTENTIALS



// should be updated to call a different function.. Coul/param.. hf etc.
double coul(int i, int Z){
// for pure coulomb.. - infinitesimal nucleus
	double vpot;
	if (i<=0){
		vpot=-Z/(0.01*r0);			// ?? or let it be infinite??
	}
	else if ((i>0)and(i<NGP)){
		vpot = -Z/r(i);
	}
	else{
		printf("FAILURE: Wrong gridpoint used for coul(i=%i)\n",i);
		return 0;
	}
	return vpot;
}


// Potential for a uniformly charged sphere
double sphere(int i, int Z, int A, double cR0){
	// cR0 is optional input of R_0 radius!
	// cR0 is input in fm (fempto m!) ~1.1A^1/3 etc etc..
	// if cR0=0, then uses default

	double vpot;
	
	//workout default value for nuclear rms radius
	double rms;
	if (A>9){
		rms=(0.836*pow(A,0.333)+0.570)*(1e-15/aB);
	}
	else{
		rms=(1.15*pow(A,0.333))*(1e-15/aB);		// no idea if this is correct
	}
	
	// If input is zero, use default value
	if (cR0==0){
		cR0=sqrt(5/3)*rms;
	}
	else {
		cR0=cR0*(1e-15/aB); 			// converts input from fm to a.u.
	}


	if (i<=0){
		vpot=-Z/(0.01*r0);			// ?? or let it be infinite??
	}
	else if ((i>0)and(i<NGP)){
		if (r(i)<=cR0){
			vpot=-(Z/(2*cR0))*(3-pow(r(i),2)/pow(cR0,2));
		}
		else {
			vpot = -Z/r(i);
		}
	}
	else{
		printf("FAILURE: Wrong gridpoint used for sphere(i=%i)\n",i);
		return 0;
	}
	return vpot;
}




// GREEN parametric potential (includes Coulomb already.. can be modified for fns)
double green(int i, int Z){
// for pure coulomb.. - infinitesimal nucleus
	double vgreen;
	double iH=4.4691;		//params for Cs
	double id=0.8967;
//	double iH=3.4811;		//params for Rb
//	double id=0.7855;
	if (i==0){
		vgreen=-Z/(0.01*r0);			// ?? or let it be infinite??
	}
	else if ((i>0)and(i<NGP)){
		vgreen = -(1/r(i))*(1+((Z-1)/(iH*(exp(r(i)/id)-1)+1)));
	}
	else{
		printf("FAILURE: Wrong gridpoint used for green(i=%i)\n",i);
		return 0;
	}
	return vgreen;
}


// TIETZ parametric potential (includes Coulomb already.. can be modified for fns)
double tietz(int i, int Z){
// for pure coulomb.. - infinitesimal nucleus
	double vtietz;
	double ig=0.2445;		//params for Cs
	double it=2.0453;
	if (i==0){
		vtietz=-Z/(0.01*r0);			// ?? or let it be infinite??
	}
	else if ((i>0)and(i<NGP)){
		vtietz = -(1/r(i))*(1+((Z-1)*exp(-ig*r(i)))/pow((1+it*r(i)),2));
	}
	else{
		printf("FAILURE: Wrong gridpoint used for tietz(i=%i)\n",i);
		return 0;
	}
	return vtietz;
}


// "Approximate" (expansion) solution for Thomas-Fermi potential for neutral atom
double TFapprox(int i, int Z){
// for pure coulomb.. - infinitesimal nucleus (can change)
// Agrees very will with GREEN parametric potential [Checked for Cs only]
// Potential corresponds to V^N-1 potential
	double x;
	double phi;
	double vTF;
	if (i==0){
//		vTF=-Z/(r0);			// ?? or let it be infinite??.. This shouldn't be used
		vTF=-Z/(0.01*r0);			// ?? or let it be infinite??.. This shouldn't be used
	}
	else if ((i>0)and(i<NGP)){
//		x=1.12950781018323*r(i)*pow(Z,1/3.);
		x=1.12950781018323*r(i)*cbrt(Z);		
		phi=1/(1+0.02747*pow(x,0.5)+1.243*pow(x,1)-0.1486*pow(x,1.5)+0.2302*pow(x,2)+0.007298*pow(x,2.5)+0.006944*pow(x,3));
//		vTF = -((Z-1)*phi+1)/r(i);
		vTF = -(0.95)*((Z-1)*phi+1)/r(i);		
	}
	else{
		printf("FAILURE: Wrong gridpoint used for TFapprox(i=%i)\n",i);
		return 0;
	}
	return vTF;
}



//*********************************************************
//*********************************************************
//*********************************************************
//			INTEGRATE
// Integrates input function f(i) wrt r, from l to m. This program contains wronskian drdt
//Uses an nquad-point quadrature formula.. for nquad=1->14 (any integer)
double integrate(double *f, int l, int m){

//Note: normally, NGP is so large, that there is no difference between nquad=1 to 14!


//re-arange limits to make less fiddly
// so I don't need to ..
	int negInt=1;
	if (l>m){
		int temp=l;
		l=m;
		m=temp;
		negInt=-1;
	}
	if (m>=NGP){m=NGP-1;} 
	if (l<0){l=0;}


// Defines the `nquad'-point quadrature integration coeficents. 
// nquad can take any integer from 1 to 14 (=1 implies trapazoid rule)
	double ic[nquad];
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
		printf("FAILURE: Wrong order for integration.. check nquad\n");
		return 0;
	}

	double dd[14]={2,2,24,24,1440,1440,120960,120960,7257600,
					7257600,958003200,958003200,5230697472000,5230697472000};


	// wronskian
	for (int i=l; i<=m; i++){
		f[i]=f[i]*drdt(i);
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
	
	// rectangle method (for tests)
	double b=0;
	for (int i=l; i<m; i++){
		b=b+f[i];
	}
	b=b*h;

	return a;
}		// END integrate




//*********************************************************
//*********************************************************
//*********************************************************
// 			Numerical DE solve by integration (adams) coeficients


//******************************************************************
// coeficients for the ADAMS-MOULTON routine
int AMcoefs(double *mia, double &mid, double &miaa){

// coefs for Adams..

	if (AMO==8){
		double tia[8]={-33953,312874,-1291214,3146338,-5033120,5595358,-4604594,4467094};
		mid=3628800;
		miaa=1070017;
		for (int i=0; i<AMO; i++){
			mia[i]=tia[i];
		}
	}
	else if(AMO==7){
		double tia[7]={1375,-11351,41499,-88547,123133,-121797,139849};
		mid=120960;
		miaa=36799;
		for (int i=0; i<AMO; i++){
			mia[i]=tia[i];
		}
	}
	else if(AMO==6){
		double tia[6]={-863,6312,-20211,37504,-46461,65112};
		mid=60480;
		miaa=19087;
		for (int i=0; i<AMO; i++){
			mia[i]=tia[i];
		}
	}
	else if(AMO==5){
		double tia[5]={27,-173,482,-798,1427};
		mid=1440;
		miaa=475;
		for (int i=0; i<AMO; i++){
			mia[i]=tia[i];
		}
	}
	else if(AMO==4){
		double tia[4]={-19,106,-264,646};
		mid=720;
		miaa=251;
		for (int i=0; i<AMO; i++){
			mia[i]=tia[i];
		}
	}
	else if(AMO==3){
		double tia[3]={1,-5,19};
		mid=24;
		miaa=9;
		for (int i=0; i<AMO; i++){
			mia[i]=tia[i];
		}
	}
	else if(AMO==2){
		double tia[2]={-1,8};
		mid=12;
		miaa=5;
		for (int i=0; i<AMO; i++){
			mia[i]=tia[i];
		}
	}
	else if(AMO==1){
		double tia[1]={1};
		mid=2;
		miaa=1;
		for (int i=0; i<AMO; i++){
			mia[i]=tia[i];
		}
	}
	else{
		printf("FAILURE: No Adams-Moulton coeficients. Check AMO\n");
		return 1;
	}
	
	return 0;
}	// END AMcoefs




//******************************************************************
// Calculates the derivative of a function!
int diff(double *f, double *deriv){

// Problem: This calcs the derivative for a few (4) points AFTER pinf; where it oscilates.
// Means the first 0 of dP won't correspond to first max of P!

// This hasn't been fully tested!

/*
df/dr = df/dt * dt/dr = (df/dt) / (dr/dt) = (df/di) * (di/dt) / (dr/dt) = 
= (df/di)  / (h * dr/dt)
coeficients from: http://en.wikipedia.org/wiki/Finite_difference_coefficient
*/


	deriv[0]= (f[1]-f[0])/(h*drdt(0));
	deriv[NGP-1]= (f[NGP-1]-f[NGP-2])/(h*drdt(NGP-1));
	
	deriv[1]= (f[2]-f[0])/(2*h*drdt(1));
	deriv[NGP-2]= (f[NGP-1]-f[NGP-3])/(2*h*drdt(NGP-2));
	
	deriv[2]= (f[0]-8*f[1]+8*f[3]-f[4])/(12*h*drdt(2));
	deriv[NGP-3]= (f[NGP-5]-8*f[NGP-4]+8*f[NGP-2]-f[NGP-1])/(12*h*drdt(NGP-3));

	deriv[3]= (-1*f[0]+9*f[1]-45*f[2]+45*f[4]-9*f[5]+1*f[6])/(60*h*drdt(3));
	deriv[NGP-4]= (-1*f[NGP-7]+9*f[NGP-6]-45*f[NGP-5]+45*f[NGP-3]-9*f[NGP-2]+1*f[NGP-1])/(60*h*drdt(NGP-4));

	for (int i=4; i<(NGP-4);i++){
//		deriv[i] = ( (1/56)*f[i-4] - (4/21)*f[i-3] + 1*f[i-2] - 4*f[i-1] 
//		- (1/56)*f[i+4] + (4/21)*f[i+3] - 1*f[i+2] + 4*f[i+1] ) / (5*h*drdt(i));
		deriv[i] = ( (1/8)*f[i-4] - (4/3)*f[i-3] + 7*f[i-2] - 28*f[i-1] 
		- (1/8)*f[i+4] + (4/3)*f[i+3] - 7*f[i+2] + 28*f[i+1] ) / (35*h*drdt(i));
	}



}


//******************************************************************
// coeficients for the OUTINT routine
int OIcoefs(double (*oie)[amo2], double *oia, double &oid){

// coefs for outint..

	if (AMO==8){
		double tie[8][8]={	-1338,2940,-2940,2450,-1470,588,-140,15,
							-240,-798,1680,-1050,560,-210,48,-5,
							60,-420,-378,1050,-420,140,-30,3,
							-32,168,-672,0,672,-168,32,-3,
							30,-140,420,-1050,378,420,-60,5,
							-48,210,-560,1050,-1680,798,240,-15,
							140,-588,1470,-2450,2940,-2940,1338,105,
							-960,3920,-9408,14700,-15680,11760,-6720,2283};
		double tia[8]={-105,15,-5,3,-3,5,-15,105};
		oid=840;
		for (int i=0; i<AMO; i++){
			oia[i]=tia[i];
			for (int j=0; j<AMO; j++){
				oie[i][j]=tie[i][j];
			}
		}
	}
	else if (AMO==7){
		double tie[7][7]={	-609,1260,-1050,700,-315,84,-10,
							-140,-329,700,-350,140,-35,4,
							42,-252,-105,420,-126,28,-3,
							-28,126,-420,105,252,-42,4,
							35,-140,350,-700,329,140,-10,
							-84,315,-700,1050,-1260,609,60,
							490,-1764,3675,-4900,4410,-2940,1089};
		double tia[7]={-60,10,-4,3,-4,10,-60};
		oid=420;
		for (int i=0; i<AMO; i++){
			oia[i]=tia[i];
			for (int j=0; j<AMO; j++){
				oie[i][j]=tie[i][j];
			}
		}
	}
	else if (AMO==6){
		double tie[6][6]={	-77,150,-100,50,-15,2,
							-24,-35,80,-30,8,-1,
							9,-45,0,45,-9,1,
							-8,30,-80,35,24,-2,
							15,-50,100,-150,77,10,
							-72,225,-400,450,-360,147};
		double tia[6]={-10,2,-1,1,-2,10};
		oid=60;
		for (int i=0; i<AMO; i++){
			oia[i]=tia[i];
			for (int j=0; j<AMO; j++){
				oie[i][j]=tie[i][j];
			}
		}
	}
	else if (AMO==5){
		double tie[5][5]={	-65,120,-60,20,-3,
							-30,-20,60,-15,2,
							15,-60,20,30,-3,
							-20,60,-120,65,12,
							75,-200,300,-300,137};
		double tia[5]={-12,3,-2,3,-12};
		oid=60;
		for (int i=0; i<AMO; i++){
			oia[i]=tia[i];
			for (int j=0; j<AMO; j++){
				oie[i][j]=tie[i][j];
			}
		}
	}
	else{
		printf("FAILURE: No Adams-Moulton (OUTINT) coeficients. Check AMO\n");
		return 1;
	}

	return 0;
}	// 	END OIcoefs







//*********************************************************
//*********************************************************
//*********************************************************
// 			OTHER


// ``Exact'' Dirac energies (for comparison) [for pure coulomb/infinitesimal nucleus]
double diracen(int z, int n, int k){
// 
	double g=sqrt(k*k-aa2*z*z);
	double diracE=c2/sqrt(1+(aa2*z*z)/pow((g+n-fabs(k)),2))-c2;	
	return diracE;
}



//******************************************************************
// Function that uses LAPACK [dgetrf+dgetri] to invert a matrix
// Optionally prints out before/after matrices, and their product (which should be 1)
extern "C" void dgetrf_(int*, int*, double*, int*, int*, int*);
extern "C" void dgetri_(int*, double*, int*, int*, double*, int*, int*);
// above definitions must be made!
int invertmat(double (*matrix)[amo2], double (*inverse)[amo2], int dim){

	int info;				 	// LAPACK output: =0 means OK
	double work;
	int ipvt[dim];				// row swaps?
	int LW=3*dim;				//dimension of `workspace'

	double tempmat[dim][dim];	

	int printmatrices=0;	//=1 to print the matrices!
	
	for (int i=0; i<dim; i++){
		for (int j=0; j<dim; j++){
			tempmat[i][j]=matrix[i][j];
		}
	}

	dgetrf_(&dim,&dim,*tempmat,&dim,ipvt,&info);
	dgetri_(&dim,*tempmat,&dim,ipvt,&work,&LW,&info); 
	
	for (int i=0; i<dim; i++){
		for (int j=0; j<dim; j++){
			inverse[i][j]=tempmat[i][j];
		}
	}
	
	// For testing the matrix inversion!
	// Prints out incoming matrix, its inverse, and their product (which should be I)
	if (printmatrices==1){
		double tempel;
		printf("\n***********Matrices********\n");
		printf("\nOriginal Matrix (A): \n");
		for (int i=0; i<dim; i++){
			for (int j=0; j<dim; j++){
				printf("% 8.1f ",matrix[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		printf("Inverted Matrix B=Inv(A): \n");
		for (int i=0; i<dim; i++){
			for (int j=0; j<dim; j++){
				printf("% 3.1e ",inverse[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		double AB[dim][dim];
		printf("Product AxB (Should be Identity): \n");
		for (int i=0; i<dim; i++){
			for (int j=0; j<dim; j++){
				tempel=0;
				for (int x=0; x<dim; x++){
					tempel=tempel+matrix[i][x]*inverse[x][j];
					AB[i][j]=tempel;
				}
				printf("% 7.1g",AB[i][j]);
			}
			printf("\n");	
		}
		printf("\n");
	} //END print matrices

	if(info!=0){
		printf("Problem with LAPACK inverting matrix! Error code: %i\n",info);	
	}		

	return info;
}	// END invertmat









