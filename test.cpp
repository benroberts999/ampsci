#include "ElectronOrbitals.h"
#include "atomInfo.h"
#include "INT_quadratureIntegration.h"
#include <iostream>


int main(void){

  clock_t ti,tf;
  ti = clock();

  printf("\nRunning SolveDBS for Local potential\n");
  printf("*******************************\n");
  //printf("NGP=%i, Size of box: Rmax=%.1f a.u., h=%f\n\n",NGP,r(NGP-1),h);

  int Z=1;
  int A=1;
  int n_max = 5;
  ElectronOrbitals wf(Z,A,2000,1.);

  wf.localBoundState(n_max);

  printf("Grid: pts=%i h=%7.5f Rmax=%5.1f\n\n",wf.ngp,wf.h,wf.r[wf.ngp-1]);


  printf(" n l_j    k  R_inf its eps   En (au)            Error (au)\n");
  int num_states = wf.nlist.size();
  for(int i=0; i<num_states; i++){
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double del = wf.en[i] - wf.diracen(wf.Z,n,k);
    double rinf = wf.r[wf.pinflist[i]];
    //std::cout<<wf.diracen(wf.Z,n,k)<<"\n";
    printf("%2i %s_%i/2 (%2i)  %3.0f %i %5.0e  %.15f  %8.1e\n",
        n,atinfo_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        wf.en[i],del);
  }

  std::cout<<"\n\n";

  // Calculate the expectation value of r^rpow for each state in list:
  int rpow=1;
  printf("\nExpectation value of r^%i\n",rpow);
  for (int s=0; s<num_states; s++){
    std::vector<double> rad;
    for (int i=0; i<wf.ngp; i++){
      //rad[i]=(P[s][i]*P[s][i]+Q[s][i]*Q[s][i])*(pow(r(i),rpow));
      double x=(wf.p[s][i]*wf.p[s][i]+wf.q[s][i]*wf.q[s][i])*pow(wf.r[i],rpow);
      rad.push_back(x);
    }
    double R=INT_integrate(rad,wf.drdt,wf.h);
    printf("<%i % i|r^%i|%i % i> = % .10f\n",wf.nlist[s],wf.klist[s],rpow,
        wf.nlist[s],wf.klist[s],R);
  }

  // HERE!? Is there some problem with the formula?? a missing n or kappa??
  // Testing Dirac Eq. by evaluating <a|H|a> - ME of Hamiltonian
  printf("\nTesting Dirac Eq: <n|H|n> (test of numerical uncertainty)\n");
  int ngp=wf.ngp;
  double alpha = wf.alpha;
  double a2 = pow(alpha,2);
  for (int s=0; s<num_states; s++){
    //std::vector<double> Ps,Qs;
    //Ps.push_back(wf.p[s]);
    //Qs.push_back(wf.p[s]);
    std::vector<double> dP(ngp),dQ(ngp);
    INT_diff(wf.q[s],dQ,wf.drdt,wf.h);
    INT_diff(wf.p[s],dP,wf.drdt,wf.h);
    for (int i=0; i<NGP; i++){
      rad[i]=(
            2*wf.p[s][i]*dQ[i]/alpha
            // P[s][i]*dQ[i]/aa
            //-Q[s][i]*dP[i]/aa
          +  ((-2*wf.klist[s])/(r[i]*alpha))*wf.p[s][i]*wf.p[s][i]
            +wf.v[i]*(wf.p[s][i]*wf.p[s][i]+wf.q[s][i]*wf.q[s][i])
           -(2./a2)*wf.q[s][i]*wf.q[s][i]
          );
    }
    double R=INT_integrate(rad,wf.drdt,wf.h);
    double fracdiff=(R-wf.en[s])/wf.en[s];
    printf("<%i% i|H|%i% i> = % .15f, E(%i% i) = % .15f; % .2e\n",wf.nlist[s],wf.klist[s],R,wf.nlist[s],wf.klist[s],wf.en[s],fracdiff);
  }




  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\n t=%.3f ms.\n",total_time);


}
