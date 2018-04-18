#include "ElectronOrbitals.h"
#include "ATI_atomInfo.h"
#include "INT_quadratureIntegration.h"
#include <iostream>
#include <fstream>

int main(void){

  clock_t ti,tf;
  ti = clock();

  int Z,A;
  int n_max,l_max;
  int ngp = 2000;
  double r0,rmax;
  double varalpha;
  int iextra;
  std::ifstream ifile;
  ifile.open("h-like.in");
  {
    std::string junk;
    ifile >> Z >> A;            getline(ifile,junk);
    ifile >> n_max >> l_max;    getline(ifile,junk);
    ifile >> r0 >> rmax >> ngp; getline(ifile,junk);
    ifile >> varalpha;          getline(ifile,junk);
    ifile >> iextra;            getline(ifile,junk);
  }
  ifile.close();
  bool extra=false;
  if(iextra==1) extra=true;

  printf("\nRunning SolveDBS for Local H-like potential, Z=%i\n",Z);
  printf("*************************************************\n");

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  if(A!=0) wf.sphericalNucleus();

  //Solve the Dirac equation for H-like ions:
  wf.hydrogenLike(n_max,l_max);

  printf("Grid: pts=%i h=%7.5f Rmax=%5.1f\n",wf.ngp,wf.h,wf.r[wf.ngp-1]);
  if(varalpha!=1) printf("varalpha = c/c_eff = %.1e  ",varalpha);
  if(varalpha<1) std::cout<<"(non-relativistic scenario)\n";
  if(varalpha>1) std::cout<<"(hyper-relativistic scenario)\n";

  std::cout<<"\n";

  printf(" n l_j    k  R_inf its eps     En (au)            Error (au)\n");
  int num_states = wf.nlist.size();
  for(int i=0; i<num_states; i++){
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double del = wf.en[i] - wf.diracen(wf.Z,n,k);
    double rinf = wf.r[wf.pinflist[i]];
    printf("%2i %s_%i/2 (%2i)  %3.0f %3i  %5.0e  %.15f  %7.0e\n",
        n,ATI_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        wf.en[i],del);
  }

  if(extra){
    // Calculate the expectation value of r^rpow for each state in list:
    printf("\nExpectation value of r^n (radial integral)\n");
    std::cout<<"          ";
    for(int in=-2; in<=2; in ++){
      if(in==0) continue;
      printf(" <nk|r^%i|nk>   ",in);
    }
    std::cout<<"\n";
    for (int s=0; s<num_states; s++){
      printf("%2i%s_%i/2 : ",wf.nlist[s],ATI_l(ATI_l_k(wf.klist[s])).c_str(),
        ATI_twoj_k(wf.klist[s]));
      for(int in=-2; in<=2; in ++){
        if(in==0) continue;
        std::vector<double> rad1;
        for (int i=0; i<wf.ngp; i++){
          double x1=(wf.p[s][i]*wf.p[s][i]+wf.q[s][i]*wf.q[s][i])*pow(wf.r[i],in);
          rad1.push_back(x1);
        }
        double R1=INT_integrate(rad1,wf.drdt,wf.h);
        printf("%13.8f, ",R1);
      }
      std::cout<<"\n";
    }

    // Testing Dirac Eq. by evaluating <a|H|a> - ME of Hamiltonian
    printf("\nTesting wavefunctions: <n|H|n>  (numerical error)\n");
    double alpha = wf.alpha;
    double a2 = pow(alpha,2);
    for (int s=0; s<num_states; s++){
      std::vector<double> dQ(wf.ngp);
      INT_diff(wf.q[s],wf.drdt,wf.h,dQ);
      std::vector<double> rad;
      for (int i=0; i<wf.ngp; i++){
        double x1=2*wf.p[s][i]*dQ[i]/alpha;
        double x2=-2*wf.klist[s]*wf.p[s][i]*wf.q[s][i]/(wf.r[i]*alpha);
        double x3=-2*pow(wf.q[s][i],2)/a2;
        double x4=wf.vnuc[i]*(pow(wf.p[s][i],2)+pow(wf.q[s][i],2));
        rad.push_back(x1+x3+x2+x4);
      }
      double R=INT_integrate(rad,wf.drdt,wf.h);
      double fracdiff=(R-wf.en[s])/wf.en[s];
      printf("<%i% i|H|%i% i> = % .15f, E(%i% i) = % .15f; % .0e\n",wf.nlist[s],
          wf.klist[s],wf.nlist[s],wf.klist[s],R,wf.nlist[s],wf.klist[s],
          wf.en[s],fracdiff);
    }
  }

  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.3f ms.\n",total_time);

  return 0;
}
