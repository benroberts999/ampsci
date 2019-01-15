#include "HartreeFockClass.h"
#include "ElectronOrbitals.h"
#include "ATI_atomInfo.h"
#include "PRM_parametricPotentials.h"
#include "WIG_369j.h"
#include "INT_quadratureIntegration.h"
#include <vector>
#include <cmath>
/*
Calculates self-consistent Hartree-Fock potential, including exchange.
Solves all core and valence states.

//NB: To-do notes int header

*/

//******************************************************************************
HartreeFock::HartreeFock(ElectronOrbitals &wf, const std::string &in_core,
  double eps_HF)
{

  m_eps_HF = eps_HF;
  if(eps_HF>1) m_eps_HF = pow(10,-1*eps_HF);//can give as log..

  //Store pointer to ElectronOrbitals object internally
  p_wf = &wf;
  m_ngp = p_wf->ngp;

  starting_approx_core(in_core);
  m_num_core_states = p_wf->num_core_states;

  //store l, 2j, and "kappa_index" in arrays for faster/easier access
  for(int kappa : p_wf->kappa){
    twoj_list.push_back(ATI::twoj_k(kappa));
    kappa_index_list.push_back(index_from_kappa(kappa));
  }

  //Run HF for all core states
  hartree_fock_core();
}


//******************************************************************************
void HartreeFock::hartree_fock_core(){

  double eta1=0.35;
  double eta2=0.7; //this value after 4 its

  form_core_Lambda_abk(p_wf->kappa);

  vex.resize(m_num_core_states,std::vector<double>(m_ngp));//Also
  initialise_m_arr_v_abk_r_core();

  std::vector<double> vdir_old(m_ngp);
  std::vector< std::vector<double> > vex_old(m_num_core_states,
    std::vector<double>(m_ngp));

  //Start the HF itterative procedure:
  int hits;
  double eta = eta1;
  for(hits=1; hits<MAX_HART_ITS; hits++){
    if(hits==4) eta = eta2;
    if(hits==32) eta = eta1; //?

    //Form new v_dir and v_ex:
    vdir_old = p_wf->vdir;
    vex_old = vex;

    form_vabk_core();
    form_vdir(p_wf->vdir,false);
    form_approx_vex_core(vex);

    for(int j=0; j<m_ngp; j++){
      p_wf->vdir[j] = eta*p_wf->vdir[j] + (1.-eta)*vdir_old[j];
      for(auto i=0ul; i<m_num_core_states; i++){
        vex[i][j] = eta*vex[i][j] + (1.-eta)*vex_old[i][j];
      }
    }

    //Solve Dirac Eq. for each state in core, using Vdir+Vex:
    std::vector<double> en_old = p_wf->en;
    double t_eps = 0;
    for(auto i=0ul; i<m_num_core_states; i++){
      //calculate de from PT
      double del_e = 0;
      for(int j=0; j<p_wf->pinflist[i]; j+=5){
        double dv = p_wf->vdir[j]+vex[i][j]-vdir_old[j]-vex_old[i][j];
        del_e += dv*(pow(p_wf->f[i][j],2) + pow(p_wf->g[i][j],2))*p_wf->drdt[j];
      }
      del_e*=p_wf->h*5;
      double en_guess = en_old[i] + del_e;
      if(en_guess>0) en_guess = en_old[i]; //safety, should never happen
      p_wf->reSolveDirac(i,en_guess,vex[i],1); //only go to 1/10^2 here
      //t_eps: weighted average of (de)/e for each orbital:
      double sfac = 2.*p_wf->kappa[i]*p_wf->occ_frac[i]; //|2k|=2j+1
      t_eps += fabs(sfac*(p_wf->en[i]-en_old[i])/en_old[i]);
    }//core states
    t_eps /= (int(p_wf->num_core_electrons)*eta);

    //Force all core orbitals to be orthogonal to each other
    p_wf->orthonormaliseOrbitals(1);

    printf("\rHF core        it:%3i eps=%6.1e              ",hits,t_eps);
    std::cout<<std::flush;
    if(t_eps<m_eps_HF && hits>2) break;
  }//hits
  std::cout<<"\n";

  //Now, re-solve core orbitals with higher precission
  for(auto i=0ul; i<m_num_core_states; i++)
    p_wf -> reSolveDirac(i,p_wf->en[i],vex[i],15);
  p_wf -> orthonormaliseOrbitals(2);
  // + re-solve direct potential (higher precission)?
  // form_vdir(p_wf->vdir,false);

}

//******************************************************************************
void HartreeFock::solveValence(int n, int kappa)
{

  twoj_list.push_back(ATI::twoj_k(kappa));
  kappa_index_list.push_back(index_from_kappa(kappa));
  extend_Lambda_abk(kappa);
  extend_m_arr_v_abk_r_valence(kappa);

  //just use direct to solve initial
  p_wf->solveLocalDirac(n,kappa,0,1);
  auto a = p_wf->nlist.size() - 1; //index of this valence state
  int twoJplus1 = ATI::twoj_k(kappa)+1; //av. over REL configs
  p_wf->occ_frac.push_back(1./twoJplus1);

  double eta1=0.35;
  double eta2=0.7; //this value after 4 its
  double eta = eta1;
  std::vector<double> vexa(m_ngp);
  std::vector<double> vexa_old(m_ngp);
  int hits;
  for(hits=1; hits<MAX_HART_ITS; hits++){
    if(hits==4) eta=eta2;
    if(hits==30) eta=eta1;
    // vexa_old = vexa;
    double en_old = p_wf->en[a];
    //Form new exch. potential:
    form_vabk_valence(a);
    vexa_old = vexa;
    form_approx_vex_a(a,vexa);
    for(int i=0; i<m_ngp; i++) vexa[i] = eta*vexa[i]+(1.-eta)*vexa_old[i];
    //Use P.T. to calculate energy change:
    double en_new = 0;
    for(int i=0; i<p_wf->pinflist[a]; i+=5) en_new += (vexa[i]-vexa_old[i])
      *(pow(p_wf->f[a][i],2) + pow(p_wf->g[a][i],2))*p_wf->drdt[i];
    en_new = p_wf->en[a] + en_new*p_wf->h*5;
    //Solve Dirac using new potential:
    p_wf->reSolveDirac(a,en_new,vexa,1);
    double eps = fabs((p_wf->en[a]-en_old)/(eta*en_old));
    //Force valence states to be orthogonal to each other + to core:
    p_wf->orthonormaliseValence(a,1);
    printf("\rHF val:%2lu %2i %2i | %3i eps=%6.1e  en=%11.8f                  "
     ,a,n,kappa,hits,eps,p_wf->en[a]);
    std::cout<<std::flush;
    if(eps<m_eps_HF && hits>2) break;
  }
  std::cout<<"\n";

  //Re-solve w/ higher precission
  p_wf->reSolveDirac(a,p_wf->en[a],vexa,15);
  p_wf->orthonormaliseValence(a,2);

}


//******************************************************************************
double HartreeFock::calculateCoreEnergy()
/*
Calculates the total HF core energy:
  E = \sum_a [ja]e_a - 0.5 \sum_(ab) (R^0_abab - \sum_k L^k_ab R^k_abba)
where:
  R^k_abcd = Integral [f_a*f_c + g_a*g_c] * v^k_bd
  R^0_abab is not absymmetric
  R^k_abba _is_ ab symmetric
*/
{
  double Etot=0;
  for(auto a=0ul; a<m_num_core_states; a++){
    double E1=0, E2=0, E3=0;
    double xtjap1 = (twoj_list[a]+1)*p_wf->occ_frac[a];
    E1 += xtjap1*p_wf->en[a];
    for(auto b=0ul; b<m_num_core_states; b++){
      double xtjbp1 = (twoj_list[b]+1)*p_wf->occ_frac[b];
      std::vector<double> &v0bb = get_v_aa0(b);
      double R0f2 = INT::integrate4(p_wf->f[a],p_wf->f[a],v0bb,p_wf->drdt);
      double R0g2 = INT::integrate4(p_wf->g[a],p_wf->g[a],v0bb,p_wf->drdt);
      E2 += xtjap1*xtjbp1*(R0f2+R0g2);
      //take advantage of symmetry for third term:
      if(b>a) continue;
      double y = (a==b) ? 1 : 2;
      int kmin = abs(twoj_list[a] - twoj_list[b])/2;
      int kmax = (twoj_list[a] + twoj_list[b])/2;
      std::vector<std::vector<double> > &vabk = get_v_abk(a,b);
      for(int k=kmin; k<=kmax; k++){
        double L_abk = get_Lambda_abk(a,b,k);
        if(L_abk==0) continue;
        int ik = k - kmin;
        double R0f3 =
          INT::integrate4(p_wf->f[a],p_wf->f[b],vabk[ik],p_wf->drdt);
        double R0g3 =
          INT::integrate4(p_wf->g[a],p_wf->g[b],vabk[ik],p_wf->drdt);
        E3 += y*xtjap1*xtjbp1*L_abk*(R0f3+R0g3);
      }
    }
    {
      Etot += E1-0.5*(E2-E3)*p_wf->h; //update running total
    }
  }
  return Etot;
}

//******************************************************************************
void HartreeFock::starting_approx_core(const std::string &in_core)
/*
Starting approx for HF. Uses Green parametric
Later, can put other options if you want.
*/
{
  p_wf -> vdir.resize(m_ngp); //make sure correct size
  int Z = p_wf -> Znuc();
  //Get default values for Green potential
  double Gh,Gd;  //Green potential parameters
  PRM::defaultGreenCore(Z,Gh,Gd);
  //Fill the the potential, using Greens PRM
  for(int i=0; i<m_ngp; i++)
    p_wf -> vdir[i] = PRM::green(Z,p_wf->r[i],Gh,Gd);
  //First step: Solve each core state using above parameteric potential
  p_wf -> solveInitialCore(in_core,1);//1 in 10
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
This routine re-sizes the m_arr_Lambda_nmk array
New routine for valence? Or make so can re-call this one??
*/
{
  m_arr_Lambda_nmk.clear(); //should already be empty!

  //Find largest existing kappa index
  int max_kappa_index = 0;
  for(auto i=0ul; i<kappa.size(); i++){
    int kappa_index = index_from_kappa(kappa[i]);
    if(kappa_index > max_kappa_index) max_kappa_index = kappa_index;
  }
  m_max_kappa_index_so_far = max_kappa_index;

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
      }//k
      Lmk.push_back(Lk);
    }//m
    m_arr_Lambda_nmk.push_back(Lmk);
  }//n

}

//******************************************************************************
void HartreeFock::extend_Lambda_abk(int kappa_a)
/*

Note: there is code overlap with: form_core_Lambda_abk
Could create a new function? Or just leave it?
Or, merge this with above (if statements etc??)
Note: don't just add this one, because might have skipped indexes!
Add all we might need, keep order matchine index!
*/
{
  int n_a = index_from_kappa(kappa_a);
  if(n_a <= m_max_kappa_index_so_far) return; //already done

  for(int n=m_max_kappa_index_so_far+1; n<=n_a; n++){
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
    m_arr_Lambda_nmk.push_back(Lmk);
  }
  m_max_kappa_index_so_far = n_a;
}

//******************************************************************************
double HartreeFock::get_Lambda_abk(unsigned long a, unsigned long b, int k) const
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

  if(n>m) return m_arr_Lambda_nmk[n][m][k-kmin];
  return m_arr_Lambda_nmk[m][n][k-kmin];
}


//******************************************************************************
void HartreeFock::initialise_m_arr_v_abk_r_core()
/*
Initialise (re-size) array to store CORE HF screening functions, v^k_ab(r)
Note: only for core. These are stored in m_arr_v_abk_r array (class member)
*/
{
  m_arr_v_abk_r.clear();
  m_arr_v_abk_r.resize(m_num_core_states);
  for(auto a=0ul; a<m_num_core_states; a++){
    m_arr_v_abk_r[a].resize(a+1);
    int tja = twoj_list[a];
    for(auto b=0ul; b<=a; b++){
      int tjb = twoj_list[b];
      int num_k = (tja>tjb) ? (tjb+1) : (tja+1);
      m_arr_v_abk_r[a][b].resize(num_k);
      for(int ik=0; ik<num_k; ik++){
        m_arr_v_abk_r[a][b][ik].resize(m_ngp);
      }//k
    }//b
  }//a
}
//******************************************************************************
void HartreeFock::extend_m_arr_v_abk_r_valence(int kappa_a)
/*
This enlargens the m_arr_v_abk_r to make room for the valence states
*/
{
  std::vector<std::vector<std::vector<double> > > v_abk_tmp(m_num_core_states);
  int tja = 2*abs(kappa_a)-1;  // |2k|=2j+1
  for(auto b=0ul; b<m_num_core_states; b++){
    int tjb = twoj_list[b];
    int num_k = (tja>tjb) ? (tjb+1) : (tja+1);
    v_abk_tmp[b].resize(num_k, std::vector<double>(m_ngp));
  }//b
  m_arr_v_abk_r.push_back(v_abk_tmp);
}

//******************************************************************************
void HartreeFock::calculate_v_abk(int a, int b, int k,
  std::vector<double> & vabk)
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
  int irmax = std::min(p_wf->pinflist[a],p_wf->pinflist[b]); //pass in instead?

  double Ax=0, Bx=0;
  for(int i=0; i<irmax; i++) Bx += p_wf->drdt[i]*
    (p_wf->f[a][i]*p_wf->f[b][i] + p_wf->g[a][i]*p_wf->g[b][i])
      /pow(p_wf->r[i],k+1);

  if(a==b) irmax = m_ngp; //For direct part, can't cut!
  vabk[0] = Bx*p_wf->h;
  for(int i=1; i<irmax; i++){
    double Fdr = p_wf->drdt[i-1]*
      (p_wf->f[a][i-1]*p_wf->f[b][i-1]+p_wf->g[a][i-1]*p_wf->g[b][i-1]);
    Ax = Ax + Fdr*pow(p_wf->r[i-1],k);
    Bx = Bx - Fdr/pow(p_wf->r[i-1],k+1);
    vabk[i] = p_wf->h*(Ax/pow(p_wf->r[i],k+1) + Bx*pow(p_wf->r[i],k));
  }
  for(int i=irmax; i<m_ngp; i++) vabk[i] = 0;
}

//******************************************************************************
void HartreeFock::form_vbb0()
/*
When doing Hartree (no exchange) only need v^0_bb
Don't call this as well as form_vabk_core, not needed (won't break though)
*/
{
  for(auto b=0u; b<p_wf->num_core_states; b++)
    calculate_v_abk(b,b,0,m_arr_v_abk_r[b][b][0]);
}

//******************************************************************************
void HartreeFock::form_vabk_core()
/*
Calculates [calls calculate_v_abk] and stores the v^k_ab Hartree-Fock sreening
 functions for (a,b) in the core.
Takes advantage of a/b symmetry. Skips if Lambda=0 (integral=0 from angles)
Note: only for core-core states! (for now?)
*/
{
  #pragma omp parallel for
  for(auto a=0ul; a<m_num_core_states; a++){
    for(auto b=0ul; b<=a; b++){
      int kmin = abs(twoj_list[a] - twoj_list[b])/2;
      int kmax = (twoj_list[a] + twoj_list[b])/2;
      for(int k=kmin; k<=kmax; k++){
        if(get_Lambda_abk(a,b,k)==0) continue;
        calculate_v_abk(a,b,k,m_arr_v_abk_r[a][b][k-kmin]);
      }//k
    }//b
  }//a
}

//******************************************************************************
void HartreeFock::form_vabk_valence(unsigned long w)
/*
Calculates [calls calculate_v_abk] and stores the Hartree-Fock screening
functions v^k_wb for a single (given) valence state (w=valence, b=core).
Stores in m_arr_v_abk_r
*/
{
  #pragma omp parallel for
  for(auto b=0ul; b<m_num_core_states; b++){
    int kmin = abs(twoj_list[w] - twoj_list[b])/2;
    int kmax = (twoj_list[w] + twoj_list[b])/2;
    for(int k=kmin; k<=kmax; k++){
      if(get_Lambda_abk(w,b,k)==0) continue;
      calculate_v_abk(w,b,k,m_arr_v_abk_r[w][b][k-kmin]);
    }//k
  }//b
}

//******************************************************************************
std::vector<std::vector<double> >& HartreeFock::get_v_abk(int a, int b)
/*
Returns a reference to a 2D-array (a subset of the m_arr_v_abk_r array)
Returned array is of form: array[ik][r]; ik runs from 0 -> |kmax-kmin+1|
  array.size() = |kmax-kmin+1| = number of k's
  array[0].size() = ngp
Allows to call for any a,b, even though only calculated for a>=b (symmetry)
*/
{
  if(a>b) return m_arr_v_abk_r[a][b];
  return m_arr_v_abk_r[b][a];
}
//******************************************************************************
std::vector<double>& HartreeFock::get_v_aa0(int a)
// Same as above, but for v^0_aa, only need to return 1D array: array[r]
// array.size()=ngp
{
  return m_arr_v_abk_r[a][a][0];
}

//******************************************************************************
void HartreeFock::form_vdir(std::vector<double> &vdir, bool re_scale)
/*
Forms the direct part of the potential.
Must call either form_vbb0 or form_vabk_core first!
Doesn't calculate, assumes m_arr_v_abk_r array exists + is up-to-date
If re_scale==true, will scale by (N-1)/N. This then given the averaged Hartree
 potential (local, same each state, no exchange). re_scale=false by default
*/
{
  for(int i=0; i<m_ngp; i++) vdir[i] = 0;
  double sf = re_scale? (1. - (1.)/p_wf->num_core_electrons) : 1;
  for(auto b=0ul; b<m_num_core_states; b++){
    double f = (twoj_list[b]+1)*p_wf->occ_frac[b];
    std::vector<double> &v0bb = get_v_aa0(b);
    for(int i=0; i<m_ngp; i++) vdir[i] += f*v0bb[i]*sf;
  }//b
}

//******************************************************************************
void HartreeFock::form_approx_vex_core(std::vector<std::vector<double> > &vex)
/*
Forms the 2D "approximate" exchange potential for each core state, a.
NOTE: Must call form_vabk_core first!
Doesn't calculate, assumes m_arr_v_abk_r array exists + is up-to-date
*/
{
  #pragma omp parallel for
  for(auto a=0ul; a<m_num_core_states; a++){
    form_approx_vex_a(a,vex[a]);
  }
}

//******************************************************************************
void HartreeFock::form_approx_vex_a(unsigned long a, std::vector<double> &vex_a)
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

  //define references to orbitals (this is just to save typing..)
  auto &fa =  p_wf->f[a];
  auto &ga =  p_wf->g[a];

  for(auto b=0ul; b<m_num_core_states; b++){ // b!=b
    if(b==a) continue;
    auto &fb =  p_wf->f[b];
    auto &gb =  p_wf->g[b];
    double x_tjbp1 = (twoj_list[b]+1)*p_wf->occ_frac[b];
    int irmax = std::min(p_wf->pinflist[a],p_wf->pinflist[b]);
    int kmin = abs(twoj_list[a] - twoj_list[b])/2;
    int kmax = (twoj_list[a] + twoj_list[b])/2;
    std::vector<std::vector<double> > &vabk = get_v_abk(a,b);
    std::vector<double> v_Fab(m_ngp);//hold "fraction" psi_a*psi_b/(psi_a^2)
    for(int i=0; i<irmax; i++){
      // This is the approximte part! Divides by psi_a
      if(fabs(fa[i])<1.e-3) continue;
      double fac_top = fa[i]*fb[i] + ga[i]*gb[i];
      double fac_bot = fa[i]*fa[i] + ga[i]*ga[i];
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
    double x_tjap1 = (twoj_list[a]+1)*p_wf->occ_frac[a];
    int irmax = p_wf->pinflist[a];
    int kmax = twoj_list[a];
    std::vector<std::vector<double> > &vaak = get_v_abk(a,a);
    for(int k=0; k<=kmax; k++){
      double L_abk = get_Lambda_abk(a,a,k);
      if(L_abk==0) continue;
      for(int i=0; i<irmax; i++) vex_a[i] += -1*L_abk*vaak[k][i]*x_tjap1;
    }//k
  }//if a in core

}
