#include "HartreeFockClass.h"
#include "ElectronOrbitals.h"
#include "ATI_atomInfo.h"
#include "PRM_parametricPotentials.h"
#include "WIG_369j.h"
#include <vector>
#include <cmath>

// XXX Rename to Coulomb Integrals??

//******************************************************************************
HartreeFock::HartreeFock(ElectronOrbitals &wf)
{

  startingApprox();

  //store l, 2j, and "kappa_index" in arrays
  for(int i=0; i<wf.num_core_states; i++){
    l_list.push_back(ATI::l_k(wf.kappa[i]));
    twoj_list.push_back(ATI::twoj_k(wf.kappa[i]));
    kappa_index_list.push_back(index_from_kappa(wf.kappa[i]));
  }

  //initialise the non-square, non-regular array:
  arr_v_abk_r.clear();
  arr_v_abk_r.resize(wf.num_core_states);
  for(int a=0; a<wf.num_core_states; a++){
    arr_v_abk_r[a].resize(a+1);
    int tja = twoj_list[a];
    for(int b=0; b<=a; b++){
      int tjb = twoj_list[b];
      int num_k = (tja>tjb) ? (tjb+1)/2 : (tja+1)/2;
      arr_v_abk_r[a][b].resize(num_k); //right? or +1? XXX
      for(int ik=0; ik<num_k; ik++){ //every second!
        arr_v_abk_r[a][b][ik].resize(wf.ngp);
      }
    }
  }

}

//******************************************************************************
//XXX Explain index here
int HartreeFock::kappa_from_index(int i){
  // (-1,n+1) * (int(n/2)+1)
  //this gives kappa, for n=0,1,2,3...
  int sign = (i%2) ? 1 : -1;//if i is even, sign is negative
  return sign*(int(i/2)+1);
}
int HartreeFock::index_from_kappa(int ka){
  if(ka>0) return 2*ka-1;
  return 2*abs(ka)-2;
}
int HartreeFock::twoj_from_index(int i){
  return 2*int(i/2)+1;
}

//******************************************************************************
void HartreeFock::form_Lambda_abk()
{

  int max_kappa_index = 0;
  for(int i=0; i<wf.num_core_states; i++){
    kappa_index = index_from_kappa(wf.kappa[i]);
    if(kappa_index > max_kappa_index) max_kappa_index = kappa_index;
  }

  arr_Lambda_nmk.clear(); //should already be empty!

  for(int n=0; n<=max_kappa_index; n++){
    int tja = twoj_from_index(n);
    std::vector<std::vector<double> > Lmk;
    for(int m=0; m<=n; m++){
      int tjb = twoj_from_index(m);
      int kmin = (tja - tjb)/2; //don't need abs, as m<=n => ja>=jb
      int kmax = (tja + tjb)/2;
      int ikmax = (kmax-kmin)/2;
      std::vector<double> Lk(ikmax+1);
      for(int ik=0; ik<=ikmax; ik++){
        int k = 2*ik + kmin;
        double tjs = WIG::threej_2(tja,tjb,2*k,-1,1,0);
        //WIG::parity(la,lb,k); //not needed! every 2nd k is good! XXX CHECK?
        Lk[ik] = tjs*tjs;
      }
      Lmk.push_back(Lk);
    }
    arr_Lambda_nmk.push_back(Lmk);
  }

}

//******************************************************************************
double Lambda_abk(int a, int b, int k) const
//gets Lambda
//a nd b and STATE indexes
//NOTE: There are a few if's here. Maybe better not?
//Write another UNSAFE verion! Can only be called w/ correct a,b,k's!! ??
{

  int n = kappa_index_list[a];
  int m = kappa_index_list[b];

  int kmin = abs(twoj_list[a] - twoj_list[b])/2;
  if(k<kmin) return 0;
  if(k%2 != kmin%2) return 0; //if kmin is even, k must be even

  int kmax = (twoj_list[a] + twoj_list[b])/2;
  if(k>kmax) return 0;

  int ik = (k-kmin)/2
  if(n<=m) return arr_Lambda_nmk[n][m][ik];
  return arr_Lambda_nmk[m][n][ik];

}
//******************************************************************************
double unsafe_Lambda_abk(int a, int b, int ik) const
//No checks at all...
//If not any faster, just kill this version!!
{
  int n = kappa_index_list[a];
  int m = kappa_index_list[b];
  //int ik = (k-kmin)/2
  if(n<=m) return arr_Lambda_nmk[n][m][ik];
  return arr_Lambda_nmk[m][n][ik];
}

//******************************************************************************
void HartreeFock::form_v_abk()
{

  for(int a=0; a<wf.num_core_states; a++){
    for(int b=0; b<=a; b++){
      int irmax = std::min(wf.pinflist[a],wf.pinflist[b]);


      int kmin = abs(twoj_list[a] - twoj_list[b])/2;
      int kmax = (twoj_list[a] + twoj_list[b])/2;
      for(int k=kmin; k<=kmax; k+=2){ //every second!

        int ik = (k-kmin)/2
        double Labk = unsafe_Lambda_abk(a,b,ik);

        double Ax=0, Bx=0;
        for(int i=0; i<irmax; i++) Bx += wf.drdt[i]*
          (wf.f[a][i]*wf.f[b][i]+wf.g[a][i]*wf.g[b][i])/pow(wf.r[i],k+1);

          arr_v_abk_r[a][b][ik][0] = 0;
          for(int i=1; i<irmax; i++){
            double Fdr = wf.drdt[i-1]*
              (wf.f[a][i-1]*wf.f[b][i-1]+wf.g[a][i-1]*wf.g[b][i-1]);
              Ax = Ax + Fdr*pow(wf.r[i-1],k);
              Bx = Bx - Fdr/pow(wf.r[i-1],k+1);
              arr_v_abk_r[a][b][ik][i]
                = wf.h*(Ax/pow(wf.r[i],k+1) + Bx*pow(wf.r[i],k));
          }
          for(int i=irmax; i<wf.ngp; i++) arr_v_abk_r[a][b][ik][i] = 0;

      }
    }
  }
}



//******************************************************************************
void HartreeFock::startingApprox()
/*
Starting approx for HF. Uses Green parametric
Later, can put other options if you want.
*/
{
  wf.vdir.resize(wf.ngp); //make sure correct size

  //Get default values for Green potential
  double Gh,Gd;  //Green potential parameters
  PRM::defaultGreenCore(wf.Z,Gh,Gd);
  //Fill the the potential, using Greens PRM
  for(int i=0; i<ngp; i++) wf.vdir[i] = PRM::green(wf.Z,wf.r[i],Gh,Gd);

  //First step: Solve each core state using above parameteric potential
  wf.solveInitialCore(2); //2, since don't need high accuray here [1 in 10^2]

}
