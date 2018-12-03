#include "HartreeFockClass.h"
#include "ElectronOrbitals.h"
#include "ATI_atomInfo.h"
#include "PRM_parametricPotentials.h"
#include "WIG_369j.h"
#include <vector>
#include <cmath>
//Add change update git

// XXX Rename to Coulomb Integrals??

//******************************************************************************
HartreeFock::HartreeFock(ElectronOrbitals &wf)
{

  //wf = in_wf;

  startingApprox(wf);
  m_ngp = wf.ngp;
  m_num_core_states = wf.num_core_states;



  //store l, 2j, and "kappa_index" in arrays
  for(int i=0; i<wf.num_core_states; i++){
    twoj_list.push_back(ATI::twoj_k(wf.kappa[i]));
    // l_list.push_back(ATI::l_k(wf.kappa[i]));
    kappa_index_list.push_back(index_from_kappa(wf.kappa[i]));
  }

  // form_Lambda_abk(wf.kappa);

  hartree_fock_core(wf);

}

//******************************************************************************
void HartreeFock::initialise_arr_v_bb0_r(){
  arr_v_bb0_r.clear();
  arr_v_bb0_r.resize(m_num_core_states, std::vector<double>(m_ngp));
}
//******************************************************************************
void HartreeFock::initialise_arr_v_abk_r(){
  arr_v_abk_r.clear();
  arr_v_abk_r.resize(m_num_core_states);
  for(int a=0; a<m_num_core_states; a++){
    arr_v_abk_r[a].resize(a+1);
    int tja = twoj_list[a];
    for(int b=0; b<a; b++){
      int tjb = twoj_list[b];
      int num_k = (tja>tjb) ? (tjb+1) : (tja+1);
      arr_v_abk_r[a][b].resize(num_k); //right? or +1? XXX
      for(int ik=0; ik<num_k; ik++){ //every second!
        arr_v_abk_r[a][b][ik].resize(m_ngp);
      }
    }
  }
}

//******************************************************************************
void HartreeFock::hartree_fock_core(ElectronOrbitals &wf){

  int Ncs = wf.num_core_states;
  double eta1=0.35;
  double eta2=0.7; //this value after 4 its
  int ngp = wf.ngp;

  form_Lambda_abk(wf.kappa);

  vex.resize(Ncs,std::vector<double>(ngp));//Also
  initialise_arr_v_bb0_r();
  initialise_arr_v_abk_r();

  //these: move to different function!!
  // form_vabk(wf); //note: call other one if HART only!

  // form_vdir(wf.vdir,wf,false);
  // form_approx_vex_core(vex,wf);

  int MAX_HART_ITS = 64;
  double eps_HF = 1.e-6; //XXX

  //Start the HF itterative procedure:
  int hits;
  double eta = eta1;

  for(hits=1; hits<MAX_HART_ITS; hits++){
    if(hits==4) eta = eta2;
    //if(hits==32) eta = eta1; //?

    //Form new v_dir and v_ex:

    std::vector<double> vdir_old = wf.vdir;
    std::vector< std::vector<double> > vex_old = vex;
    form_vabk(wf);
    //form_vbb0(wf);
    form_vdir(wf.vdir,wf,false);
    form_approx_vex_core(vex,wf);

    // std::cout<<"\n";
    // for(int j=0; j<wf.ngp; j++){
    //   std::cout<<wf.r[j]<<" "<<vdir_old[j]<<" "<<wf.vdir[j]<<"\n";
    // }

    // form_approx_vex_core(vex,wf);
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
      for(int j=0; j<wf.ngp; j++){
        double dv = wf.vdir[j]+vex[i][j]-vdir_old[j]-vex_old[i][j];
        del_e += dv*(pow(wf.f[i][j],2) + pow(wf.g[i][j],2))*wf.drdt[j];
      }
      // for(int j=0; j<wf.ngp; j++){
      //   std::cout<<wf.r[j]<<" "<<wf.vdir[j]<<" "<<53./wf.r[j]<<" "<<vex[i][j]<<"\n";
      // }
      // break;
      del_e*=wf.h;
      std::cout<<"\n"<<en_old[i]<<" "<<del_e<<" ";
      double en_guess = en_old[i] + del_e;
      if(en_guess>0) en_guess = en_old[i]; //safety, should never happen
      wf.reSolveDirac(i,en_guess,vex[i],1); //only go to 1/10^2 here
      //wf.reSolveDirac(i,en_guess,3); //only go to 1/10^2 here
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
    //std::cout<<wf.en[0]<<" ";
    std::cin.get();
  }
  std::cout<<"\n";

  //Now, re-solve core orbitals with higher precission
  // + re-solve direct potential (higher precission)
  for(int i=0; i<Ncs; i++) wf.reSolveDirac(i,wf.en[i],vex[i],15);
  wf.orthonormaliseOrbitals(1);
  // form_vdir(wf.vdir,wf,false); //+ new v_ex??
  form_vabk(wf);
  form_vdir(wf.vdir,wf,false);
  // form_approx_vex_core(vex,wf);
  for(int i=0; i<Ncs; i++) wf.reSolveDirac(i,wf.en[i],vex[i],15);
  wf.orthonormaliseOrbitals(2);
}
//******************************************************************************
void HartreeFock::startingApprox(ElectronOrbitals &wf)
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
  wf.solveInitialCore(1); //2, since don't need high accuray here [1 in 10^2]

}

//******************************************************************************
//XXX Explain index here
int HartreeFock::kappa_from_index(int i){
  // (-1,n+1) * (int(n/2)+1)
  //this gives kappa, for n=0,1,2,3...
  int sign = (i%2) ? 1 : -1;//if i is even, sign is negative
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
void HartreeFock::form_Lambda_abk(const std::vector<int> &kappa)
// XXX For now, just core.
// New routine for valence? Or make so can re-call this one??
// THIS routine re-sizes the arr_Lambda_nmk array
{

  int max_kappa_index = 0;
  for(size_t i=0; i<kappa.size(); i++){
    int kappa_index = index_from_kappa(kappa[i]);
    if(kappa_index > max_kappa_index) max_kappa_index = kappa_index;
  }
  m_max_kappa_index_core = max_kappa_index;

  //xxx careful. For valence, this needs updating!
  arr_Lambda_nmk.clear(); //should already be empty!

  for(int n=0; n<=max_kappa_index; n++){
    int tja = twoj_from_index(n);
    int la  = l_from_index(n);
    std::vector<std::vector<double> > Lmk;
    for(int m=0; m<=n; m++){
      int tjb = twoj_from_index(m);
      int lb  = l_from_index(m);
      int kmin = (tja - tjb)/2; //don't need abs, as m<=n => ja>=jb
      int kmax = (tja + tjb)/2;
      // std::cout<<"\n j,kmin: "<<0.5*tja<<" "<<0.5*tjb<<" : "<<kmin<<" -> "<<kmax<<"\n";
      int num_k = (kmax-kmin)+1; ///2;
      //std::cout<<"\n\n"<<kmin<<" "<<kmax<<" \n";
      std::vector<double> Lk(num_k); //+1 ?? XXX
      for(int k=kmin; k<=kmax; k++){
        int ik = k-kmin;
        int p = WIG::parity(la,lb,k);
        // std::cout<<"\n\n"<<kappa_from_index(n)<<" "<<kappa_from_index(m)<<" "<<k<<" "<<
        //     WIG::parity(ATI::l_k(kappa_from_index(n)),ATI::l_k(kappa_from_index(m)),k)<<"\n\n";
        // if(p==0) continue;
        double tjs = WIG::threej_2(tja,tjb,2*k,-1,1,0);
        Lk[ik] = tjs*tjs*p;
      }
      // for(int ik=0; ik<=ikmax; ik++){
      //   int k = 2*ik + kmin;
      //   double tjs = WIG::threej_2(tja,tjb,2*k,-1,1,0);
      //   //WIG::parity(la,lb,k); //not needed! every 2nd k is good! XXX CHECK?
      //   std::cout<<"\n\n"<<kappa_from_index(n)<<" "<<kappa_from_index(m)<<" "<<k<<" "<<
      //     WIG::parity(ATI::l_k(kappa_from_index(n)),ATI::l_k(kappa_from_index(m)),k)<<"\n\n";
      //   Lk[ik] = tjs*tjs;
      // }
      Lmk.push_back(Lk);
    }
    arr_Lambda_nmk.push_back(Lmk);
  }

}

// //******************************************************************************
// int HartreeFock::get_num_ks(int a, int b) const
// //NOTE: can only be called AFTER arr_v_abk_r is formed!
// {
//   if(a==b) return 1;
//   if(a>b) return arr_v_abk_r[a][b].size();
//   return arr_v_abk_r[b][a].size();
// }

//******************************************************************************
double HartreeFock::get_Lambda_abk(int a, int b, int k) const
//gets Lambda
//a nd b and STATE indexes
//NOTE: There are a few if's here. Maybe better not?
{

  int n = kappa_index_list[a];
  int m = kappa_index_list[b];

  int kmin = abs(twoj_list[a] - twoj_list[b])/2;
  if(k<kmin) return 0;
  //if(k%2 != kmin%2) return 0; //if kmin is even, k must be even

  int kmax = (twoj_list[a] + twoj_list[b])/2;
  if(k>kmax) return 0;

  int ik = (k-kmin);
  if(n>=m) return arr_Lambda_nmk[n][m][ik];
  return arr_Lambda_nmk[m][n][ik];

}

//******************************************************************************
void HartreeFock::calculate_v_abk(const ElectronOrbitals &wf,
  int a, int b, int k, std::vector<double> & vabk)
/*
Calculalates v^k_ab screening function.
Note: only call for b<=a, and for k's with non-zero angular coefs!
Since v_ab = v_ba !

Stores in vabk - must already be sized corectly!

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

  int irmax = std::min(wf.pinflist[a],wf.pinflist[b]); //? called too often?
  //could do k loop in here too? Nah
  //irmax = m_ngp; //XXX XXX XXX XXX

  double Ax=0, Bx=0;
  for(int i=0; i<irmax; i++) Bx += wf.drdt[i]*
    (wf.f[a][i]*wf.f[b][i]+wf.g[a][i]*wf.g[b][i])/pow(wf.r[i],k+1);

  if(a==b) irmax = m_ngp; //XXX
  vabk[0] = Bx;
  for(int i=1; i<irmax; i++){
    double Fdr = wf.drdt[i-1]*
      (wf.f[a][i-1]*wf.f[b][i-1]+wf.g[a][i-1]*wf.g[b][i-1]);
    Ax = Ax + Fdr*pow(wf.r[i-1],k);
    Bx = Bx - Fdr/pow(wf.r[i-1],k+1);
    vabk[i] = wf.h*(Ax/pow(wf.r[i],k+1) + Bx*pow(wf.r[i],k));
  }
  for(int i=irmax; i<wf.ngp; i++) vabk[i] = 0;

}

//******************************************************************************
void HartreeFock::form_vbb0(const ElectronOrbitals &wf)
// ONLY call this one if not doing HF!
{
  // arr_v_bb0_r //must already be correct size!
  for(int b=0; b<wf.num_core_states; b++){
    calculate_v_abk(wf,b,b,0,arr_v_bb0_r[b]);
  }
}

//******************************************************************************
void HartreeFock::form_vabk(const ElectronOrbitals &wf)
//Does k=0 case too, so don't call form_vbb0 as well!
{
  // arr_v_abk_r //must already be correct size!
  #pragma omp parallel for
  for(int a=0; a<m_num_core_states; a++){
    for(int b=0; b<a; b++){
      int kmin = abs(twoj_list[a] - twoj_list[b])/2;
      int kmax = (twoj_list[a] + twoj_list[b])/2;
      for(int k=kmin; k<=kmax; k++){
        int ik = (k-kmin);
        //XXX put check: don't calc if Lam = 0!!
        if(get_Lambda_abk(a,b,k)==0) continue;
        calculate_v_abk(wf,a,b,k,arr_v_abk_r[a][b][ik]);
      }
    }
    calculate_v_abk(wf,a,a,0,arr_v_bb0_r[a]);
  }
}

//******************************************************************************
std::vector<double>& HartreeFock::get_v_abk(int a, int b, int k)
//XXX make const?
//XXX Put 'pragma' debug in here - turned on/off @ compile! (no slow-down!)
//XXX Check if k out of bounds??
{
  if(a==b) return arr_v_bb0_r[a];

  //I might already know kmin! XXX
  int kmin = abs(twoj_list[a] - twoj_list[b])/2;
  int ik = (k-kmin);

  if(a>b) return arr_v_abk_r[a][b][ik];
  return arr_v_abk_r[b][a][ik]; //check! Is a copy hapenning?? XXX

  // if(k<kmin) return 0;
  // if(k%2 != kmin%2) return 0; //if kmin is even, k must be even
  // int kmax = (twoj_list[a] + twoj_list[b])/2;
  // if(k>kmax) return 0;

}


//******************************************************************************
void HartreeFock::form_vdir(std::vector<double> &vdir,
  const ElectronOrbitals &wf, bool re_scale)
{

  //clear old values
  for(int i=0; i<m_ngp; i++) vdir[i] = 0;

  double sf = re_scale? (1. - (1.)/wf.num_core_electrons) : 1;

  //put OMP here, with critical?
  //OR omp below?
  //OR diff loop structure??
  for(int b=0; b<m_num_core_states; b++){
    //int rmax = wf.pinflist[b];
    double f = (twoj_list[b]+1)*wf.occ_frac[b];
    std::vector<double> &v0bb = get_v_abk(b,b,0); // is there a copy?????
    for(int i=0; i<m_ngp; i++){
      vdir[i] += f*v0bb[i]*sf;
    }
  }

}

//******************************************************************************
void HartreeFock::form_approx_vex_core(std::vector<std::vector<double> > &vex,
  const ElectronOrbitals &wf)
{
  #pragma omp parallel for
  for(int a=0; a<m_num_core_states; a++){
    form_approx_vex_a(a,vex[a],wf); //XXX remember to size first!
  }
}

//******************************************************************************
void HartreeFock::form_approx_vex_a(int a, std::vector<double> &vex_a,
  const ElectronOrbitals &wf)
{

  for(int i=0; i<m_ngp; i++) vex_a[i] = 0;

  //Can play w/ OMP + critical here..or not (shouldn't matter!)
  for(int b=0; b<m_num_core_states; b++){
    if(b==a) continue;
    double x_tjbp1 = (twoj_list[b]+1)*wf.occ_frac[b]; //what about ocf[a]??? XXX
    int irmax = std::min(wf.pinflist[a],wf.pinflist[b]);
    //irmax = m_ngp; // XXX XXX XXX NGP
    //int Nk = get_num_ks(a,b);
    int kmin = abs(twoj_list[a] - twoj_list[b])/2;
    int kmax = (twoj_list[a] + twoj_list[b])/2;
    for(int k=kmin; k<=kmax; k++){
      //XXX some stuff in here that doesn't depend on k!! take outside of loop??
      double f1 = get_Lambda_abk(a,b,k);
      if(f1==0) continue;
      double f2 = -1*(x_tjbp1)*f1;
      std::vector<double> &vabk = get_v_abk(a,b,k); // is there a copy?????
      //std::cout<<"\n\n"<<b<<" "<<k<<" "<<f1<<"\n\n";
      for(int i=0; i<irmax; i++){
        if(fabs(wf.f[a][i])<1.e-3) continue;
        //double Fab=0;
          double fac_top = wf.f[a][i]*wf.f[b][i] + wf.g[a][i]*wf.g[b][i];
          double fac_bot = wf.f[a][i]*wf.f[a][i] + wf.g[a][i]*wf.g[a][i];
          double Fab = fac_top/fac_bot;
        vex_a[i] += f2*vabk[i]*Fab;
      }
    }
  }
  //now, do a=b, k=0 case - ONLY if a is in the core!
  if(a<m_num_core_states){
    double x_tjbp1 = (twoj_list[a]+1)*wf.occ_frac[a]; //what about ocf[a]??? XXX
    //int irmax = wf.pinflist[a];
    //irmax = m_ngp; //XXX
    std::vector<double> &vabk = get_v_abk(a,a,0); // is there a copy?????
    double f1 = -1*x_tjbp1*get_Lambda_abk(a,a,0);
    for(int i=0; i<m_ngp; i++) vex_a[i] += f1*vabk[i];
  }

}
