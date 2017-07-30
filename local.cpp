/*---------------------
Thu 22 Jan 2015 17:34:30 AEDT

Routing for Local Potential (coulomb, parametric)
Calls "solveBS"

This is the file that uses solvebs.cpp program to solve dirac equation for local potentials.

It should take input from "local.in" file


---------------------*/
//using namespace std;
#include "funs.h"
#include "solvebs.h"
#include "params.h"





//#include "./solvebs.cpp"   // this program solves dirac Eq. (Contains solveBS)



//**************************************************************************
//**************************************************************************
//**************************************************************************
//**************************************************************************




//**************************************************************************
//**************************************************************************
//**************************************************************************
//**************************************************************************
int main(void){

  clock_t t;
  t = clock();

  printf("\nRunning SolveDBS for Local potential\n");
  printf("*******************************\n");
  printf("NGP=%i, Size of box: Rmax=%.1f a.u., h=%f\n\n",NGP,r(NGP-1),h);


/*
//should be input!
  int Z=1;
  int n=1;
  int ka=-1;
  double en=-0.1;

  // if a zero (or positive) input, then Guess!
  // EDIT this, to work with parametric pot!
  if(en>=0){
    en=-0.5*(pow(Z,2)/pow(n,2));
  }

*/

  // XXX Have a look-up table, so can input at. symbol instead of Z
  // (Entering Z should be OK too!)
  // Then, if I enter A=0, it uses the default A, otherwise the given A
  // (I think A is more common than N)

  //input parameters
  int Z=0;
  int Za=0;     //Is this A?? or N??

  //read in from input file:
  std::ifstream myfile;
  myfile.open ("local.in");
        myfile >> Z;
        myfile >> Za;
  myfile.close();

  printf("%Z=i, A(?)=%i\n",Z,Za);



  int pinf,its;
  double eps;
  double en;
  int ka;

  // arrays for the wavefunction. (atm alpha inside Q)
  double p[NGP],q[NGP];


  const int MAXSTATES=50; //???

  //Work out the number of states needed:
  int nmin=1;
  int nmax=4;
  int kmax=4;    //kmax=4 means up to f-states, =3 => d states etc.  |ka+1|<=kmax
  int nlist=0;
  for (int i=nmin; i<=nmax; i++){
    if (i<=kmax){
      nlist=nlist+2*i-1;
    }
    else {
      nlist=nlist+7;
    }
  }


  // Arrays for storing the wavefunctions
  double P[MAXSTATES][NGP],Q[MAXSTATES][NGP];
  double ev[MAXSTATES];
  int kappa[MAXSTATES];
  int nn[MAXSTATES];

  //Fills in the potential (should be dependent on input)
  double v[NGP];
  for (int i=0; i<NGP; i++){
    v[i]=coul(i,Z);        //atm only use coulomb pot
//    v[i]=TFapprox(i,Z);
//    v[i]=sphere(i,Z,5,0);    //Uniformly charged sphere (i,Z,A,cR0)
  }

  printf("n  ka Rinf its    eps       en (a.u.)        Error\n");
  double gen;  //"Guess" energy (?)
  int a=0;     //(state number - just a label)
  for (int n=nmin; n<=nmax; n++){
    ka=0;
    for (int k=0; k<2*n-1; k++){  //k<2*n-1 does all.. k<7 upto f states
      if (ka==0){
        ka=-1;
      } else if (ka<0){
        ka=ka*-1;
      } else {
        ka=-1*(ka+1);
      }
      if (fabs(ka+1)>kmax) break;

      // Energy guess for H-like
      en=diracen(Z,n,ka);

      //Solve bound-state problem:
      solveDBS(p,q,v,Z,n,ka,en,pinf,its,eps);

      //store resultant wavefunctions:
      for (int i=0; i<NGP; i++){
        P[a][i]=p[i];
        Q[a][i]=q[i];
      }

      //store the parameters (kappa, energy, n)
      kappa[a]=ka;
      ev[a]=en;
      nn[a]=n;

      //fractional difference between calculated
      // and "Exact" Dirac energy.
      double fracdiff=(en-diracen(Z,n,ka))/en;

      //print results to screen: (can do this aftern?)
      printf("%2i % 3i % 4.0f % 3i % 8.1e % .15f % .2e",n,ka,r(pinf),its,eps,
             en,fracdiff);
      if (fabs(fracdiff)>1e-9){
        printf(" **\n");
      }else{
        printf("\n");
      }

      a++; //increment 'state' counter
    }//END loop over ang momentum Q#s, kappa
  }//End loop over principal Q number, n

  double rad[NGP];
  int rpow=1;    //power of r for radial integral


// Calculate the expectation value of r^rpow for each state in list:
  printf("\nExpectation value of r^%i\n",rpow);
  for (int s=0; s<nlist; s++){
    for (int i=0; i<NGP; i++){
      rad[i]=(P[s][i]*P[s][i]+Q[s][i]*Q[s][i])*(pow(r(i),rpow));
    }
    double R=integrate(rad,0,NGP-1);
    printf("<%i % i|r^%i|%i % i> = % .14f\n",nn[s],kappa[s],rpow,nn[s],kappa[s],R);
  }


  // Calculate the radial integrals of r^rpow for each state in list (no double counts):
  // still increases like 0.5*nmax^2
  printf("\nRadial integrals for r^%i\n",rpow);
  for (int s=0; s<5; s++){
    for (int b=0; b<s; b++){
      for (int i=0; i<NGP; i++){
        rad[i]=(P[s][i]*P[b][i]+Q[s][i]*Q[b][i])*(pow(r(i),rpow));
      }
      double R=integrate(rad,0,NGP-1);
      printf("<%i % i|r^%i|%i % i> = % .14f\n",nn[s],kappa[s],rpow,nn[b],kappa[b],R);
    }
  }



// HERE!? Is there some problem with the formula?? a missing n or kappa??
  // Testing Dirac Eq. by evaluating <a|H|a> - ME of Hamiltonian
  printf("\nTesting Dirac Eq: <n|H|n> (test of numerical uncertainty)\n");
  for (int s=0; s<nlist; s++){
    for (int i=0; i<NGP; i++){
      rad[i]=(
            ((-2*kappa[s])/(r(i)*aa))*P[s][i]*Q[s][i]
            +v[i]*(P[s][i]*P[s][i]+Q[s][i]*Q[s][i])
           -(2/aa2)*Q[s][i]*Q[s][i]
          );
    }
    double R=integrate(rad,0,NGP-1);
    double fracdiff=(R-ev[s])/ev[s];
    printf("<%i% i|H|%i% i> = % .15f, E(%i% i) = % .15f; % .2e\n",nn[s],kappa[s],nn[s],kappa[s],R,nn[s],kappa[s],ev[s],fracdiff);
  }




/*
  // Orthogonality tests:
  printf("\nTesting orthormality [Radial Integrals only!]\n");
  for (int s=0; s<nlist; s++){
    for (int b=0; b<nlist; b++){
      if (kappa[s]==kappa[b]) {
        for (int i=0; i<NGP; i++){
          rad[i]=(P[s][i]*P[b][i]+Q[s][i]*Q[b][i]);
        }
        double R=integrate(rad,0,NGP-1);
        printf("<%i % i|%i % i> = % .5e\n",nn[s],kappa[s],nn[b],kappa[b],R);
      }
    }
  }
*/



/*
// Test derivatives.

  double dP[NGP],ddP[NGP],dQ[NGP],PP[NGP],QQ[NGP];
  for(int i=0;i<NGP;i++){
    PP[i]=P[0][i];
    QQ[i]=Q[0][i];
  }
  diff(PP,dP);
  diff(QQ,dQ);
  diff(dP,ddP);

//  for(int i=0;i<NGP;i++){
//    printf("P[%.2e]=%.2e,  dP[%.2e]=%.2e,  ddP[%.2e]=%.2e\n",r(i),PP[i],r(i),dP[i],r(i),ddP[i]);
//  }





  // TEST Dirac equation with derivatives! This should give the energy!


//VP + c dQ/dr - c k/r Q = en P
//(V-2c^2)Q - c dP/dr - c k/r P = en Q


  double dirE;
  for (int i=1000; i<3000; i++){
    dirE=(v[i]*PP[i] + cc*dQ[i] - (cc*kappa[0]/r(i))*QQ[i])-ev[0]*PP[i];
    printf("%i, %.2e, %.10f\n",i,r(i),dirE);
  }

*/


  /*

//----- FOR Cesium
  int Z=5;
  double v[NGP];
  for (int i=0; i<NGP; i++){
//    v[i]=green(i,Z);      //atm only use coulomb pot
    v[i]=TFapprox(i,Z);
//    v[i]=green(i,Z)-coul(i,Z)+sphere(i,Z,133,0);
//    v[i]=tietz(i,Z);
//    v[i]=tietz(i,Z)-coul(i,Z)+sphere(i,Z,133,0);
//    v[i]=(1/3.)*tietz(i,Z)+(1/3.)*green(i,Z)+(1/3.)*TFapprox(i,Z);
//    v[i]=(1/3.)*tietz(i,Z)+(1/3.)*green(i,Z)+(1/3.)*TFapprox(i,Z)-coul(i,Z)+sphere(i,Z,133,0);
  }

  printf(" n  ka  Rinf its  eps     En0    en (a.u.)          level (/cm)\n");
  double gen;
  int kmax=2;    //kmax=4 means up to f-states, =3 => d states etc.  |ka+1|<=kmax
  for (int n=2; n<4; n++){
    ka=0;
    for (int k=0; k<2*n-1; k++){  //k<2*n-1 does all.. k<7 upto f states
      if (ka==0){
        ka=-1;
      } else if (ka<0){
        ka=ka*-1;
      }
      else {
        ka=-1*(ka+1);
      }
      if (fabs(ka+1)>kmax){break;}

      // Energy guess for Cesium!
      int el;
      if (ka>0){
        el=ka;
      }else{
        el=-ka-1;
      }
      int ncore=1; //=5 for Cs etc
      gen=-0.2/((n-ncore)*sqrt(el+1)); // works for Cs up to f states for n<17
      if(n>16){gen=gen/2;}
      en=gen;

      solveDBS(p,q,v,Z,n,ka,en,pinf,its,eps);
      double levelCs=(en+0.143576576173033)*invcm;
      printf("%2i % 3i % 4.0f % 3i % 8.1e % .3f % .15f % .2f\n",n,ka,r(pinf),its,eps,gen,en,levelCs);
    }
  }


  */


  t = clock() - t;
  printf ("\nThis took me %li clicks (%f seconds).\n",t,((double)t)/CLOCKS_PER_SEC);



}




















