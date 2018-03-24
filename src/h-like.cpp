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
  std::ifstream ifile;
  ifile.open("h-like.in");
  {
    std::string junk;
    ifile >> Z >> A;            getline(ifile,junk);
    ifile >> n_max >> l_max;    getline(ifile,junk);
    ifile >> r0 >> rmax >> ngp; getline(ifile,junk);
    ifile >> varalpha;          getline(ifile,junk);
  }
  ifile.close();

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


  // Calculate the expectation value of r^rpow for each state in list:
  printf("\nExpectation value of r^n\n");
  for (int s=0; s<num_states; s++){
    std::vector<double> rad1,rad2;
    for (int i=0; i<wf.ngp; i++){
      double x1=(wf.p[s][i]*wf.p[s][i]+wf.q[s][i]*wf.q[s][i])*pow(wf.r[i],1);
      double x2=(wf.p[s][i]*wf.p[s][i]+wf.q[s][i]*wf.q[s][i])*pow(wf.r[i],-1);
      rad1.push_back(x1);
      rad2.push_back(x2);
    }
    double R1=INT_integrate(rad1,wf.drdt,wf.h);
    double R2=INT_integrate(rad2,wf.drdt,wf.h);
    printf("<%i% i|r^n|%i% i> n=1 %13.10f, n=-1 %13.10f\n",wf.nlist[s],
        wf.klist[s],wf.nlist[s],wf.klist[s],R1,R2);
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

  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.3f ms.\n",total_time);

  std::cout<<"\n\n";
  std::vector<double> pc(ngp),qc(ngp);
  std::vector<double> pc2(ngp),qc2(ngp);
  double ec = 0.2;
  solveContinuum(pc,qc,ec,wf.vnuc,wf.Z,-1,wf.r,wf.drdt,wf.h,wf.ngp,wf.alpha);
  solveContinuum(pc2,qc2,0.4,wf.vnuc,wf.Z,-1,wf.r,wf.drdt,wf.h,wf.ngp,wf.alpha);

  std::vector<double> p1p2(ngp);
  for(int i=0; i<wf.ngp; i++)
    p1p2[i]=pc[i]*pc2[i] + qc[i]*qc2[i];

  double del=INT_integrate(p1p2,wf.drdt,wf.h,0,wf.ngp);
  std::cout<<del<<"\n";


  std::ofstream ofile;
  ofile.open("cont.txt");
  for(int i=0; i<wf.ngp; i++){
    ofile<<wf.r[i]<<" "<<pc[i]<<" "<<pc2[i]<<"\n";
  }
  ofile.close();

  // Find the r's for psi=zero, two consec => period
  //Once period is converged enough, can normalise by comparison with
  // exact (asymptotic) solution (??)
  // Find "maximum" amplitude, by using a quadratic fit to 2 nearest points
  // Scale by ratio of this maximum to max of analytic soln
  double xa=1,xb=pc[wf.ngp-300];
  double wk1=-1, wk2=0;
  for(int i=wf.ngp-300; i<wf.ngp; i++){
    xa=xb;
    xb=pc[i];
    if(xb*xa<0){
      double r1 = (wf.r[i]*pc[i-1]-wf.r[i-1]*pc[i])/(pc[i-1]-pc[i]);
      double ya=xb,yb=xb;
      for(int j=i+1; j<wf.ngp; j++){
        ya=yb;
        yb=pc[j];
        if(ya*yb<0){
          double r2 = (wf.r[j]*pc[j-1]-wf.r[j-1]*pc[j])/(pc[j-1]-pc[j]);
          std::cout<<i<<" "<<j<<" "<<r2-r1<<" "<<0.5*pow(3.1416/(r2-r1),2)<<"\n";
          wk1 = wk2;
          wk2 = r2-r1;
          //std::cout<<wk1<<" "<<wk2<<" "<<
          break;
        }
      }
      if(fabs(wk1-wk2)<1.e-4) break;
    }

    /*
    psi ~ A Cos(kr)
    r0 is 'zero', r1, r2 are first points on either side
    y1 = psi(r1) etc
    then:
    A = y1 - dy * dr1^2 / (dr2^2 - dr1^2)
    dy  = y2-y1
    dr1 = r1-r0
    dr2 = r2-r0
    [also have converse...if not the same, take average?]
    A = y2 - dy * dr2^2 / (dr2^2 - dr1^2)
    XXX NOTE! haven't checked this!! Check using Mathematica!! XXX
    XXX NO! doesn't seem to work... XXX
    */

  }












  return 0;
}
