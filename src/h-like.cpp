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

//ngp=500000;



//   //double out = fitQuadratic(1007.,1008.,1009.,1.39351,1.40745,1.2818,true);
//   //double out = fitQuadratic(1007.5,1007.6,1007.7,1.41817,1.41886,1.41813,true);
//   //double out = fitQuadratic(1007.0,1007.5,1008.,1.39351,1.41817,1.40745,true);
//   double out = fitQuadratic(.0,.5,.1,6.,7.,5.);
//   std::cout<<out<<"\n";
//   printf("%.8f\n",out);
//
//
// return 1;






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







  //std::cout<<"\n\n";


  double ec = 0.25;
  int kc=-1;

  double Zion=1.;
  std::vector<double> rc=wf.r;
  std::vector<double> drdtc=wf.drdt;
  int NGPc=wf.ngp;
  double last_r = wf.r[wf.ngp-1];
  std::vector<double> vc=wf.vnuc;

  double lam = 1.e7; //XXX ???
  double r_asym = (Zion + sqrt(4.*lam*ec+pow(Zion,2)))/(2.*ec);

  //XXX testing only!!!
  double hc = wf.h;
  double hc_target = (M_PI/15)/sqrt(2.*ec);
  if(hc>hc_target){
    printf("WARNING: Grid not dense enough for ec=%.3f (h=%.2f, need h<%.2f)",
      ec,hc,hc_target);
    if(hc>2*hc_target){
      printf("FAILURE: Grid not dense enough for ec=%.3f (h=%.2f, need h<%.2f)",
      ec,hc,hc_target);
      return 1;
    }
  }

  //Set up temporary continuum grid:
  // Extend grid for continuum state. Goes to ~20% more than 'anymp region'
  while(true){
    double r_new = last_r + hc;
    rc.push_back(r_new);
    drdtc.push_back(hc/wf.h); //??? correct??
    vc.push_back(-Zion/r_new);
    NGPc++;
    last_r = r_new;
    if(last_r>1.2*r_asym) break;
  }

  //Solve Dirac equation on temporary (extended) grid
  std::vector<double> pc(NGPc),qc(NGPc);
  solveContinuum(pc,qc,ec,vc,wf.Z,kc,rc,drdtc,wf.h,NGPc,wf.alpha);


  std::ofstream ofile;
  ofile.open("cont.txt");
  for(int i=0; i<NGPc; i++){
    ofile<<rc[i]<<" "<<pc[i]<<" "<<qc[i]<<"\n";
  }
  ofile.close();

  //find index of 'asymptotically large r' region
  int i0 = wf.ngp;
  for(int i=wf.ngp; i<NGPc; i++){
    if(rc[i]>r_asym){
      i0=i-1;
      break;
    }
  }
  //std::cout<<"\nr_asym="<<r_asym<<", i0="<<i0<<"\n\n";


  // Find the r's for psi=zero, two consec => period
  //Once period is converged enough, can normalise by comparison with
  // exact (asymptotic) solution (??)
  //nb: looks for convergence between r(NGP) and r_asym
  // If doesn't 'converge' in this region, uses r_asym
  double xa=1,xb=pc[wf.ngp];
  double wk1=-1, wk2=0;
  for(int i=wf.ngp; i<i0; i++){
    xa=xb;
    xb=pc[i];
    if(xb*xa<0){
      //Use linear extrapolation to find exact r of the zero:
      double r1 = (rc[i]*pc[i-1]-rc[i-1]*pc[i])/(pc[i-1]-pc[i]);
      double ya=xb,yb=xb;
      for(int j=i+1; j<NGPc; j++){
        ya=yb;
        yb=pc[j];
        if(ya*yb<0){
          //Use linear extrapolation to find exact r of the next zero:
          double r2 = (rc[j]*pc[j-1]-rc[j-1]*pc[j])/(pc[j-1]-pc[j]);
          wk1 = wk2;
          wk2 = r2-r1;
          break;
        }
      }
      if(fabs(wk1-wk2)<1.e-4){
        //check for "convergence"
        i0=i;
        break;
      }
    }
  }
  //XXX Note: this method only 'kinda' works???
  // a) am I not going out far enough?
  // b) Or, it's just not that accurate? Seems acurate to 1 - 0.1% ?


  // Find "maximum" amplitude, by using a quadratic fit to 2 nearest points
  // Scale by ratio of this maximum to max of analytic soln
  int ntry=0, maxtry=5;
  double amp=0;
  while(ntry<maxtry){
    //find first zero after r_asym
    for(int i=i0; i<NGPc; i++){
      if(pc[i]*pc[i-1]<0){
        i0=i;
        break;
      }
    }
    //find max:
    int imax=0;
    double y0,y1,y2,y3,y4;
    double x0,x1,x2,x3,x4;
    for(int i=i0+1; i<NGPc; i++){
      if(fabs(pc[i])<fabs(pc[i-1])){
        imax=i-1;
        y0=fabs(pc[i-3]);
        y1=fabs(pc[i-2]);
        y2=fabs(pc[i-1]);
        y3=fabs(pc[i]);
        y4=fabs(pc[i+1]);
        x0=rc[i-3];
        x1=rc[i-2];
        x2=rc[i-1];
        x3=rc[i];
        x4=rc[i+1];
        break;
      }
    }
    i0++;
    ntry++;
    double out1 = fitQuadratic(x1,x2,x3,y1,y2,y3);
    double out2 = fitQuadratic(x0,x2,x4,y0,y2,y4);
    amp+=0.5*(out1+out2);
  }
  amp/=maxtry;

  // Calculate normalisation coeficient, D, and re-scaling factor:
  // D = Sqrt[alpha/(pi*eps)] <-- Amplitude of large-r p(r)
  // eps = Sqrt[en/(en+2mc^2)]
  double al2 = pow(wf.alpha,2);
  double ceps = sqrt(ec/(ec*al2+2.)); // c*eps = eps/alpha
  //printf("%.10f\n",ceps);
  double D = 1./sqrt(M_PI*ceps);
  //printf("%.10f\n",D);
  double sf=D/amp; //re-scale factor
  //Normalise the wfs:
  for(int i=0; i<NGPc; i++){
    // XXX here: transfere into 'true' p,q (from temp ones!)
    pc[i] *= sf;
    qc[i] *= sf;
  }





  //std::ofstream ofile;
  ofile.open("contN.txt");
  for(int i=0; i<NGPc; i++){
    ofile<<rc[i]<<" "<<pc[i]<<" "<<qc[i]<<"\n";
  }
  ofile.close();






  return 0;
}
