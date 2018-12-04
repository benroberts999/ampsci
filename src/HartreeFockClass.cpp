#include "HartreeFockClass.h"
#include "ElectronOrbitals.h"
#include "ATI_atomInfo.h"
#include "PRM_parametricPotentials.h"
#include "WIG_369j.h"
#include <vector>
#include <cmath>
// #include <fstream>

//******************************************************************************
HartreeFock::HartreeFock(ElectronOrbitals &wf, double eps_HF)
{

  starting_approx_core(wf);
  m_ngp = wf.ngp;
  m_num_core_states = wf.num_core_states;

  //store l, 2j, and "kappa_index" in arrays for faster/easier access
  for(int i=0; i<wf.num_core_states; i++){
    twoj_list.push_back(ATI::twoj_k(wf.kappa[i]));
    kappa_index_list.push_back(index_from_kappa(wf.kappa[i]));
  }

  //Run HF for all core states
  hartree_fock_core(wf, eps_HF);

}

//******************************************************************************
void HartreeFock::initialise_m_core_v_abk_r()
/*
Initialise (re-size) array to store CORE HF screening functions, v^k_ab(r)
Note: only for core. These are stored in m_core_v_abk_r array (class member)
*/
{
  m_core_v_abk_r.clear();
  m_core_v_abk_r.resize(m_num_core_states);
  for(int a=0; a<m_num_core_states; a++){
    m_core_v_abk_r[a].resize(a+1);
    int tja = twoj_list[a];
    for(int b=0; b<=a; b++){
      int tjb = twoj_list[b];
      int num_k = (tja>tjb) ? (tjb+1) : (tja+1);
      m_core_v_abk_r[a][b].resize(num_k);
      for(int ik=0; ik<num_k; ik++){
        m_core_v_abk_r[a][b][ik].resize(m_ngp);
      }//k
    }//b
  }//a
}

//******************************************************************************
void HartreeFock::hartree_fock_core(ElectronOrbitals &wf, double eps_HF){

  int Ncs = wf.num_core_states;
  double eta1=0.35;
  double eta2=0.7; //this value after 4 its
  int ngp = wf.ngp;

  form_core_Lambda_abk(wf.kappa);

  vex.resize(Ncs,std::vector<double>(ngp));//Also
  initialise_m_core_v_abk_r();

  //Start the HF itterative procedure:
  int hits;
  double eta = eta1;

  for(hits=1; hits<MAX_HART_ITS; hits++){
    if(hits==4) eta = eta2;
    if(hits==32) eta = eta1; //?

    //Form new v_dir and v_ex:
    std::vector<double> vdir_old = wf.vdir;
    std::vector< std::vector<double> > vex_old = vex;
    form_vabk_core(wf);
    form_vdir(wf.vdir,wf,false);
    form_approx_vex_core(vex,wf);

    for(int j=0; j<wf.ngp; j++){
      wf.vdir[j] = eta*wf.vdir[j] + (1.-eta)*vdir_old[j];
      for(int i=0; i<Ncs; i++){
        vex[i][j] = eta*vex[i][j] + (1.-eta)*vex_old[i][j];
      }
    }

    //Solve Dirac Eq. for each state in core, using Vdir+Vex:
    std::vector<double> en_old = wf.en;
    double t_eps = 0;
    for(int i=0; i<Ncs; i++){
      //calculate de from PT
      double del_e = 0;
      for(int j=0; j<wf.pinflist[i]; j+=5){
        double dv = wf.vdir[j]+vex[i][j]-vdir_old[j]-vex_old[i][j];
        del_e += dv*(pow(wf.f[i][j],2) + pow(wf.g[i][j],2))*wf.drdt[j];
      }
      del_e*=wf.h*5;
      double en_guess = en_old[i] + del_e;
      if(en_guess>0) en_guess = en_old[i]; //safety, should never happen
      wf.reSolveDirac(i,en_guess,vex[i],1); //only go to 1/10^2 here
      //t_eps: weighted average of (de)/e for each orbital:
      double sfac = 2.*wf.kappa[i]*wf.occ_frac[i]; //|2k|=2j+1
      t_eps += fabs(sfac*(wf.en[i]-en_old[i])/en_old[i]);
    }
    t_eps /= (wf.num_core_electrons*eta);

    //Force all core orbitals to be orthogonal to each other
    wf.orthonormaliseOrbitals(1);

    printf("\rHF core        it:%3i eps=%6.1e              ",hits,t_eps);
    std::cout<<std::flush;
    if(t_eps<eps_HF && hits>2) break;
  }
  std::cout<<"\n";

  //Now, re-solve core orbitals with higher precission
  for(int i=0; i<Ncs; i++) wf.reSolveDirac(i,wf.en[i],vex[i],15);
  wf.orthonormaliseOrbitals(2);
  // + re-solve direct potential (higher precission)?
  // form_vdir(wf.vdir,wf,false);


}
//******************************************************************************
void HartreeFock::starting_approx_core(ElectronOrbitals &wf)
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
  for(int i=0; i<wf.ngp; i++) wf.vdir[i] = PRM::green(wf.Z,wf.r[i],Gh,Gd);

  //First step: Solve each core state using above parameteric potential
  wf.solveInitialCore(1);//1 in 10

}

//******************************************************************************
/*   Kappa Index:
For easy array access, define 1-to-1 index for each kappa:
kappa: -1  1 -2  2 -3  3 -4  4 ...
index:  0  1  2  3  4  5  6  7 ...
kappa(i) = (-1,i+1)*(int(i/2)+1)
*/
int HartreeFock::kappa_from_index(int i){
  int sign = (i%2) ? 1 : -1;
  return sign*(int(i/2)+1);
}
int HartreeFock::index_from_kappa(int ka) const{
  if(ka>0) return 2*ka-1;
  return 2*abs(ka)-2;
}
int HartreeFock::twoj_from_index(int i) const{
  return 2*int(i/2)+1;
}
int HartreeFock::l_from_index(int i) const{
  bool even = (i%2) ? false : true;
  if(even) return i/2;
  return (i+1)/2;
}

//******************************************************************************
void HartreeFock::form_core_Lambda_abk(const std::vector<int> &kappa)
/*
Calculate + store the angular coefifienct:
  Lambda^k_ab := 3js((ja,jb,k),(-1/2,1/2,0))^2*parity(la+lb+k)
Lambda^k_ab = Lambda^k_ka,kb :: i.e., only depends on kappa!
This routine re-sizes the m_core_Lambda_nmk array
New routine for valence? Or make so can re-call this one??
*/
{
  m_core_Lambda_nmk.clear(); //should already be empty!

  //Find largest existing kappa index
  int max_kappa_index = 0;
  for(size_t i=0; i<kappa.size(); i++){
    int kappa_index = index_from_kappa(kappa[i]);
    if(kappa_index > max_kappa_index) max_kappa_index = kappa_index;
  }

  for(int n=0; n<=max_kappa_index; n++){
    int tja = twoj_from_index(n);
    int la  = l_from_index(n);
    std::vector<std::vector<double> > Lmk;
    for(int m=0; m<=n; m++){
      int tjb = twoj_from_index(m);
      int lb  = l_from_index(m);
      int kmin = (tja - tjb)/2; //don't need abs, as m\leq n => ja\geq jb
      int kmax = (tja + tjb)/2;
      std::vector<double> Lk(kmax-kmin+1,0);
      for(int k=kmin; k<=kmax; k++){
        int ik = k-kmin;
        if(WIG::parity(la,lb,k)==0) continue;
        double tjs = WIG::threej_2(tja,tjb,2*k,-1,1,0);
        Lk[ik] = tjs*tjs;
      }
      Lmk.push_back(Lk);
    }
    m_core_Lambda_nmk.push_back(Lmk);
  }

}

//******************************************************************************
double HartreeFock::get_Lambda_abk(int a, int b, int k) const
/*
Simple routine to (semi-)safely return Lambda_abk
Note: input a and b are regular ElectronOrbitals state indexes
*/
{
  //Get the kappa index for each state:
  int n = kappa_index_list[a];
  int m = kappa_index_list[b];

  int kmin = abs(twoj_list[a] - twoj_list[b])/2;

  //these just safety..can remove if sped-up!
  if(k<kmin) return 0;
  int kmax = (twoj_list[a] + twoj_list[b])/2;
  if(k>kmax) return 0;

  if(n>m) return m_core_Lambda_nmk[n][m][k-kmin];
  return m_core_Lambda_nmk[m][n][k-kmin];
}

//******************************************************************************
void HartreeFock::calculate_v_abk(const ElectronOrbitals &wf,
  int a, int b, int k, std::vector<double> & vabk)
/*
Calculalates v^k_ab screening function.
Note: should only call for a>=b, and for k's with non-zero angular coefs
(nothing bad will happen othersie, but no point!)
Since v_ab = v_ba

Stores in vabk (reference to whatever) - must already be sized corectly!

r_min := min(r,r')
rho(r') := fa(r')*fb(r') + ga(r')gb(r')
v^k_ab(r) = Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
          = Int_0^r [r'^k/r^(k+1)]*rho(r') dr'
            + Int_r^inf [r^k/r'^(k+1)]*rho(r') dr'
         := A(r)/r^(k+1) + B(r)*r^k
A(r0)  = 0
B(r0)  = Int_0^inf [r^k/r'^(k+1)]*rho(r') dr'
A(r_n) = A(r_{n-1}) + (rho(r_{n-1})*r_{n-1}^k)*dr
B(r_n) = A(r_{n-1}) + (rho(r_{n-1})/r_{n-1}^(k+1))*dr
v^k_ab(rn) = A(rn)/rn^(k+1) + B(rn)*rn^k
*/
{
  int irmax = std::min(wf.pinflist[a],wf.pinflist[b]); //pass in instead?

  double Ax=0, Bx=0;
  for(int i=0; i<irmax; i++) Bx += wf.drdt[i]*
    (wf.f[a][i]*wf.f[b][i]+wf.g[a][i]*wf.g[b][i])/pow(wf.r[i],k+1);

  if(a==b) irmax = m_ngp; //For direct part, can't cut!
  vabk[0] = Bx*wf.h;
  for(int i=1; i<irmax; i++){
    double Fdr = wf.drdt[i-1]*
      (wf.f[a][i-1]*wf.f[b][i-1]+wf.g[a][i-1]*wf.g[b][i-1]);
    Ax = Ax + Fdr*pow(wf.r[i-1],k);
    Bx = Bx - Fdr/pow(wf.r[i-1],k+1);
    vabk[i] = wf.h*(Ax/pow(wf.r[i],k+1) + Bx*pow(wf.r[i],k));
  }
  for(int i=irmax; i<m_ngp; i++) vabk[i] = 0;
}

//******************************************************************************
void HartreeFock::form_vbb0(const ElectronOrbitals &wf)
/*
When doing Hartree (no exchange) only need v^0_bb
*/
{
  for(int b=0; b<wf.num_core_states; b++){
    calculate_v_abk(wf,b,b,0,m_core_v_abk_r[b][b][0]);
  }
}

//******************************************************************************
void HartreeFock::form_vabk_core(const ElectronOrbitals &wf)
/*
Calculates [calls calculate_v_abk] and stores the v^k_ab Hartree-Fock sreening
 functions for (a,b) in the core.
Takes advantage of a/b symmetry. Skips if Lambda=0 (integral=0 from angles)
Note: only for core-core states! (for now?)
*/
{
  #pragma omp parallel for
  for(int a=0; a<m_num_core_states; a++){
    for(int b=0; b<=a; b++){
      int kmin = abs(twoj_list[a] - twoj_list[b])/2;
      int kmax = (twoj_list[a] + twoj_list[b])/2;
      for(int k=kmin; k<=kmax; k++){
        if(get_Lambda_abk(a,b,k)==0) continue;
        calculate_v_abk(wf,a,b,k,m_core_v_abk_r[a][b][k-kmin]);
      }
    }
  }
}

//******************************************************************************
std::vector<std::vector<double> >& HartreeFock::get_v_abk(int a, int b)
/*
Returns a reference to a 2D-array (a subset of the m_core_v_abk_r array)
Returned array is of form: array[ik][r]; ik runs from 0 -> |kmax-kmin+1|
  array.size() = |kmax-kmin+1| = number of k's
  array[0].size() = ngp
Allows to call for any a,b, even though only calculated for a>=b (symmetry)
*/
{
  if(a>b) return m_core_v_abk_r[a][b];
  return m_core_v_abk_r[b][a];
}
//******************************************************************************
std::vector<double>& HartreeFock::get_v_aa0(int a)
// Same as above, but for v^0_aa, only need to return 1D array: array[r]
// array.size()=ngp
{
  return m_core_v_abk_r[a][a][0];
}

//******************************************************************************
void HartreeFock::form_vdir(std::vector<double> &vdir,
  const ElectronOrbitals &wf, bool re_scale)
/*
Forms the direct part of the potential.
Must call either form_vbb0 or form_vabk_core first!
Doesn't calculate, assumes m_core_v_abk_r array exists + is up-to-date
If re_scale==true, will scale by (N-1)/N. This then given the averaged Hartree
 potential (local, same each state, no exchange). re_scale=false by default
*/
{
  for(int i=0; i<m_ngp; i++) vdir[i] = 0;
  double sf = re_scale? (1. - (1.)/wf.num_core_electrons) : 1;
  for(int b=0; b<m_num_core_states; b++){
    double f = (twoj_list[b]+1)*wf.occ_frac[b];
    std::vector<double> &v0bb = get_v_aa0(b);
    for(int i=0; i<m_ngp; i++) vdir[i] += f*v0bb[i]*sf;
  }//b
}

//******************************************************************************
void HartreeFock::form_approx_vex_core(std::vector<std::vector<double> > &vex,
  const ElectronOrbitals &wf)
/*
Forms the 2D "approximate" exchange potential for each core state, a.
NOTE: Must call form_vabk_core first!
Doesn't calculate, assumes m_core_v_abk_r array exists + is up-to-date
*/
{
  #pragma omp parallel for
  for(int a=0; a<m_num_core_states; a++){
    form_approx_vex_a(a,vex[a],wf);
  }
}

//******************************************************************************
void HartreeFock::form_approx_vex_a(int a, std::vector<double> &vex_a,
  const ElectronOrbitals &wf)
/*
Forms the 2D "approximate" exchange potential for given core state, a.
Does the a=b case seperately, since it's a little simpler
Approximate:
In order to approximate solution to HF equations, I form "local" ex. potential
  [v_ex*psi_a](r) = \sum_b v_ex^(a,b)(r) * psi_b(r)
v_ex is non-local; cannot write: [v_ex*psi_a](r) =/= v_ex(r)*psi_a(r)
Instead: define local approx: vex_a
  vex_a = [v_ex*psi_a](r) *(psi_a/psi_a^2)
        = \sum_b v_ex^(a,b)(r)*psi_b(r) * (psi_a/psi_a^2)
        = \sum_b v_ex^(a,b)(r)*(psi_b(r)*psi_a) / psi_a^2
This vex_a is then a local potential (different for each state!) that can
be used as an addition to local direct potential to solve Dirac Eq. as normal.
In theory, this is exact. Clearly, however, there is an issue when psi_a is
small. Luckily, however, we don't care as much when psi_a is small! Also,
since v_ex is already small (compared to vdir), we can make good approximation.
Therefore, I only calculate vex_a when a=b, or when |psi_a| > 1.e3
Further, largest part of v_ex is when a=b. In this case, the factor=1 is exact!
*/
{
  for(int i=0; i<m_ngp; i++) vex_a[i] = 0;

  for(int b=0; b<m_num_core_states; b++){ // b!=b
    if(b==a) continue;
    double x_tjbp1 = (twoj_list[b]+1)*wf.occ_frac[b];
    int irmax = std::min(wf.pinflist[a],wf.pinflist[b]);
    int kmin = abs(twoj_list[a] - twoj_list[b])/2;
    int kmax = (twoj_list[a] + twoj_list[b])/2;
    std::vector<std::vector<double> > &vabk = get_v_abk(a,b);
    std::vector<double> v_Fab(m_ngp);//hold "fraction" psi_a*psi_b/(psi_a^2)
    for(int i=0; i<irmax; i++){
      // This is the approximte part! Divides by psi_a
      if(fabs(wf.f[a][i])<1.e-3) continue;
      double fac_top = wf.f[a][i]*wf.f[b][i] + wf.g[a][i]*wf.g[b][i];
      double fac_bot = wf.f[a][i]*wf.f[a][i] + wf.g[a][i]*wf.g[a][i];
      v_Fab[i] = -1.*x_tjbp1*fac_top/fac_bot;
    }//r
    for(int k=kmin; k<=kmax; k++){
      double L_abk = get_Lambda_abk(a,b,k);
      if(L_abk==0) continue;
      for(int i=0; i<irmax; i++){
        if(v_Fab[i]==0) continue;
        vex_a[i] += L_abk*vabk[k-kmin][i]*v_Fab[i];
      }//r
    }//k
  }//b
  //now, do a=b, ONLY if a is in the core!
  if(a<m_num_core_states){
    double x_tjap1 = (twoj_list[a]+1)*wf.occ_frac[a];
    int irmax = wf.pinflist[a];
    int kmax = twoj_list[a];
    std::vector<std::vector<double> > &vaak = get_v_abk(a,a);
    for(int k=0; k<=kmax; k++){
      double L_abk = get_Lambda_abk(a,a,k);
      if(L_abk==0) continue;
      for(int i=0; i<irmax; i++) vex_a[i] += -1*L_abk*vaak[k][i]*x_tjap1;
    }//k
  }//if a in core

}
