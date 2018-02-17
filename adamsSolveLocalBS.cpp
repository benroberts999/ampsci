/*
Thu 06 Nov 2014 23:29:57 AEDT

- IT WORKS.... (AMO=8 fails sometimes)

Program to solve single-electron bound-state Dirac problem for a (given)
local, central potential.
Based on method presented in book by W. Johnson.
Employs the Adams-Moulton method.

solveDBS is the main routine that is called from elsewhere.
All other functions called by solveDBS.
*/
#include "adamsSolveLocalBS.h"



bool dodebug=false; //if true, will print progress messages.

/*
To do :: Jan 2018
 * smarter variable names!
 * p,q -> f,g! Check the q cancellation!
 * updated ctp matching thing! [see below]

*/


//******************************************************************************
int solveDBS(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int Z, int n, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int &pinf, int &its, double &eps, double alpha)
/*
Solves local, spherical bound state dirac equation using Adams-Moulton method.
Based on method presented in book by W. R. Johnson:
  W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007)

See also:
 * https://en.wikipedia.org/wiki/Linear_multistep_method
 * Hairer, Ernst; Nørsett, Syvert Paul; Wanner, Gerhard (1993),
   Solving ordinary differential equations I: Nonstiff problems (2nd ed.),
   Berlin: Springer Verlag, ISBN 978-3-540-56670-0.
 * Quarteroni, Alfio; Sacco, Riccardo; Saleri, Fausto (2000),
   Matematica Numerica, Springer Verlag, ISBN 978-88-470-0077-3.
 * http://mathworld.wolfram.com/AdamsMethod.html

Rough description of method:
1. Start with initial 'guess' of energy XXX* (see note below)
2. Find the "practical infinity" (psi~0), and the Classical turning point [e=v]
3. Performs 'inward' integration (Adams Moulton). Integrates from the
   practical infinity inwards to the classical turning point (ctp).
4. Performs 'outward' integration (Adams Moulton). Integrates from 0
   outwards to the classical turning point.
5. Matches the two functions at ctp, by re-scaling the 'inward' solution (for P)
6. Checks the number of nodes the wf has. If too many or too few nodes, makes
   a large change to the energy and tries again (from step 2).
   If the correct number of nodes, uses perturbation theory to make minor
   corrections to the energy to 'zoom in', then re-starts from step 2.
   Continues until this energy adjustment falls below a prescribed threshold.

XXX is it possible to  slightly change the method, so that inint and outint
go _past_ ctp (by some set number of points), and then we try to match all of
these points (instead of just 1)??

XXX Go back and re-understand how the large energy changes work, AND how
the minor (P.T.) changes work!

*/
{


// XXX
// At the moment, it appears that initial energy guess in an input.
// This is good, but should also have an option that allows an initial input
// of zero, in which case this program will make the initial guess
// XXX



  // bound state wavefunctions (Adams-moul)
  const double delep=5e-15;		//PRIMARY convergence parameter for bound state energy	(10^-11)
  const double deles=1e-11;		//SECONDAY convergence parameter for bound state energy	(X)
  const int ntry=100;			// Number of failed attempts at converging (sove bs) before error quit (30)
  const double alr=800;			// ''assymptotically large r [not what this is..]''  (=800)
  const double lde=0.2;		//amount to vary energy by for 'large' variations (0.1 => 10%)
  // XXX where to put these parameters?

  //Checks to see if legal n is requested. If not, increases n, re-calls
  //Should I do this? If there's a logic problem, probably better to know?
  if ( (fabs(ka)<=n) && (ka!=n) ){
    //XXX swap around!
    if(dodebug) printf("\nRunning SolveDBS for state %i %i:\n",n,ka);
  }
  else{
    //XXX No..shouldn't do this!
    printf("\nSate %i %i does not exist..\n",n,ka);
    // n=n+1;
    // en=-0.5*(pow(Z,2)/pow(n,2));
    // solveDBS(p,q,v,Z,n,ka,en,pinf,its,eps);
    return 1;
  }


  // Find 'l' from 'kappa' (ang. momentum Q number)
  // This is used to calculate number of nodes wf should have
  int ll;   //L ang. mom QN
  if (ka>0) ll=ka;
  else ll=-ka-1;
  //Number of nodes wf should have:
  int inodes=n-ll-1;

  //Some parameters used by the Adams Moulton method:
  int more=0,less=0;          //params for checking nodes
  double eupper=0,elower=0;   //params for """ and varying energy
  double anorm=0;             // normalisation constant
  int ctp;                    //classical turning point
  bool converged=false;
  double deltaEn=0;
  its=0;            //numer of iterations (for this n,ka)

  while (!converged){

    //Find the practical infinity 'pinf'
    //Step backwards from the last point (NGP-1) until
    // (V(r) - E)*r^2 >  alr    (alr = "asymptotically large r)
    pinf=NGP-1;
    while ( (en-v[pinf])*pow(r[pinf],2) + alr < 0 ){
      pinf=pinf-1;
    }
    if(pinf==NGP-1)
      printf("WARNING: pract. inf. = size of box for %i %i\n",n,ka);
    if(dodebug)
      printf("Practical infinity (i=%i): Pinf=%.1f a.u.\n",pinf,r[pinf]);


    //Find classical turning point 'ctp'
    //Step backwards from the "practical infinity" until
    //  V(r) > E        [nb: both V and E are <0]
    ctp=pinf;
    int d_ctp=2; //XXX here. Num points ctp +/- d_ctp. Make input!
    while ( (en-v[ctp]) < 0 ){
      ctp=ctp-1;
      if (ctp-d_ctp<=0){
        //fails if ctp<0, (or ctp>pinf?)
        printf("FAILURE 96 in solveDBS: No classical region?\n");
        return 1;
      }
    }
    if(ctp+d_ctp>=pinf){
      //Didn't find ctp! Does this ever happen?
      printf("FAILURE: Turning point at or after pract. inf. \n");
      printf("ctp=%i, d_ctp=%i, pinf=%i, NGP=%i\n",ctp,d_ctp,pinf,NGP);
      // printf("ctp(%i)=%.1f, d_ctp=%i, pinf(%i)=%.1f\n",ctp,r[ctp],d_ctp,pinf,r[pinf]);
      return 1;
    }
    if(dodebug) printf("Classical turning point (i=%i): ctp=%.1f a.u.\n",
                ctp,r[ctp]);
    if(dodebug) printf("%i %i: Pinf= %.1f,  en= %f\n",n,ka,r[pinf],en);

    //Perform the "inwards integration":
    inwardAM(p,q,en,v,ka,r,drdt,h,NGP,ctp-d_ctp,pinf,alpha);

    //save the values of wf at ctp from the 'inward' ind
    //XXX ctp update! XXX
    double ptp=p[ctp];
    double qtp=q[ctp];

    //Perform the "outwards integration"
    outwardAM(p,q,en,v,Z,ka,r,drdt,h,NGP,ctp+d_ctp,alpha);

    //Scales the inward solution to match the outward solution (for P)
    double rescale=p[ctp]/ptp;
    for (int i=ctp+1; i<=pinf; i++){
      p[i]=rescale*p[i];
      q[i]=rescale*q[i];
    }
    // The qtp value is scaled too; used in perturbation correction.
    qtp=rescale*qtp;


    //Count the number of nodes (zeros) the wf has.
    //Just counts the number of times wf changes sign
    int nozeros=0;
    double sp=p[1];
    double spn;
    for (int i=2; i<pinf; i++){
      spn=p[i];
      if (sp*spn<0){
        nozeros=nozeros+1;
      }
      if(spn!=0){
        sp=spn;
      }
    }
    if(dodebug) printf("Nodes=%i\n",nozeros);


    //checks to see if there are too many/too few nodes.
    //Makes large adjustment to energy until in correct region
    //Then, once in 'correct' region, uses perturbation theory to make small
    // changes to the energy
    double etemp;
    if(nozeros>inodes){
      //Too many nodes:
      more=more+1;
      if ((more==1)||(en<eupper)){
        eupper=en;
      }
      etemp=(1+lde)*en;
      if (less!=0){
        if(etemp<elower){
          etemp=0.5*(eupper+elower);
        }
      }
      deltaEn=fabs((en-etemp)/en);
      en=etemp;
    }else if (nozeros<inodes){
      //too few nodes:
      less=less+1;
      if ((less==1)||(en>elower)){
        elower=en;
      }
      etemp=(1-lde)*en;
      if (more!=0){
        if(etemp>eupper){
          etemp=0.5*(eupper+elower);
        }
      }
      deltaEn=fabs((en-etemp)/en);
      en=etemp;
    }else{
      // correct number of nodes.
      //From here, use perturbation theory to fine-time the energy
      if(dodebug){printf("Correct number of nodes, starting P.T.\n");}
      //double ppqq[NGP];
      std::vector<double> ppqq(NGP);
      for (int i=0; i<=pinf; i++){
        ppqq[i]=p[i]*p[i]+q[i]*q[i];    // XXX add alpha here if need!
      }
      //anorm=integrate(ppqq,0,pinf);
      anorm=INT_integrate(ppqq,drdt,h,0,pinf);
      if(dodebug) printf("anrom=%.5f\n",anorm);
      double de=  (1./alpha) * p[ctp] * (qtp-q[ctp]) / anorm ;
      deltaEn=fabs(de/en);
      etemp = en + de;
      if(dodebug){
        printf("de=%.3e, en=%.5f, et=%.5f, el=%.5f, e=%.5f\n",
          de,en,etemp,elower,eupper);
      }
      if ((less!=0)&&(etemp<elower)){
        deltaEn=fabs((en-0.5*(en+elower))/en);
        en=0.5*(en+elower);
      }
      else if ((more!=0)&&(etemp>eupper)){
        deltaEn=((en-0.5*(en+eupper))/en);
        en=0.5*(en+eupper);
      }
      else if (deltaEn<delep){
        en=etemp;
        eps=deltaEn;
        converged=true;
      }
      else {
        eps=deltaEn;
        en=etemp;
      }
    }// END: if(nozeros>inodes

    its++; //increment 'number of iterations' counter

    if(dodebug) printf("Itteration number %i,  en= %f\n",its,en);
    if(its>ntry){
      if(deltaEn<deles){
        if(dodebug) printf("Wavefunction %i %i didn't fully converge after "
                    "%i iterations, but OK.\n",n,ka,ntry);
        eps=deltaEn;
        converged=true; //kind-of
      }else{
        printf("Wavefunction %i %i didn't converge after %i itterations.\n",n,
               ka,ntry);
        if(dodebug) printf("%i %i: Pinf= %.1f,  en= %f\n",n,ka,r[pinf],en);
        eps=deltaEn;
        return 2;
      }
    }

  }// END: while (not converged)


//  if(dodebug){
//    // move to class! [as an 'info' function! - used for de-bugging]
//    int iat1=log((1+r0)/r0)/h;
//    int iat10=log((10+r0)/r0)/h;
//    printf("Radial grid sepparations at 1 a.u., 10 a.u., ctp and pinf:\n");
//    printf("At    r=1    a.u.;  dr=%.4f\n",r[iat1]-r[iat1-1]);
//    printf("At    r=10   a.u.;  dr=%.4f\n",r[iat10]-r[iat10-1]);
//    printf("ctp:  r=%3.1f a.u.;  dr=%.4f\n",r[ctp],r[ctp]-r[ctp-1]);
//    printf("pinf: r=%3.1f a.u.;  dr=%.4f\n",r[pinf],r[pinf]-r[pinf-1]);
//  }


  //normalises the wavefunction
  double an= 1/sqrt(anorm);
  if(dodebug) printf("An=%.5e\n",an);
  for (int i=0; i<=pinf; i++){
    p[i]=an*p[i];
    q[i]=an*q[i];
  }
  for (int i=(pinf+1); i<NGP; i++){
    // this prob. not neccissary, but kills remainders
    p[i]=0;
    q[i]=0;
  }

  if(dodebug)
      printf("NGP=%i, Size of box: Rmax=%.1f a.u., h=%f\n",NGP,r[NGP-1],h);

  if(dodebug){
    printf("Converged for %i %i to %.3e after %i iterations;\n",n,ka,eps,its);
    printf("%i %i: Pinf(%i)=%.1f,  \nMy energy:  en= %.15f\n",n,ka,pinf,
           r[pinf],en);
    // double fracdiff=(en-diracen(Z,n,ka))/en;
    // printf("Exact Dirac:  Ed= %.15f; Del=%.2g%% \n\n",
    //     diracen(Z,n,ka),fracdiff);
  }

  return 0;
}  // END solveDBS




//*********************************************************
//*********************************************************
//*********************************************************




//******************************************************************************
int outwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int Z, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int ctp, double alpha)
/*
Program to start the OUTWARD integration.
Starts from 0, and uses an expansion(?) to go to (nol*AMO).
Then, it then call ADAMS-MOULTON, to finished (from nol*AMO+1 to ctp)

XXX rename to outwardAM ? amOutInt?

*/
{

  const int nol=1;				// # of outdir runs [finds first nol*AMO+1 points (3)]

  double az = Z*alpha;                  //Z*alpha
  double c2 = 1./pow(alpha,2);
  double ga=sqrt(pow(ka,2)-pow(az,2));  //'gamma' factor


  //initial wf values
  // P(r) = r^gamma u(r)
  // Q(r) = r^gamma v(r)
  double u0=1;
  double v0;
  if (ka>0){
    v0=-(ga+ka)/az;
  }
  else {
    v0=az/(ga-ka);
  }
  p[0]=0;
  q[0]=0;

  //Coeficients used by the method
  //double ie[amo2][amo2];
  //double ia[amo2];
  std::vector< std::vector<double> > ie(AMO,std::vector<double>(AMO));
  std::vector<double> ia(AMO);
  double id;
  getOutwardCoefs(ie,ia,id);

  // loop through and find first nol*AMO points of wf
  for (int ln=0; ln<nol; ln++){
    int i0=ln*AMO+1;

    //defines/populates coefs
    //double coefa[amo2],coefb[amo2],coefc[amo2],coefd[amo2];
    //double em[amo2][amo2];
    std::vector<double> coefa,coefb,coefc,coefd;
    std::vector< std::vector<double> > em(AMO,std::vector<double>(AMO));
    for (int i=0; i<AMO; i++){
      // coefa[i]=-id*h*(ga+ka)*dror(i+i0);
      // coefb[i]=-id*h*(en+2*c2-v[i+i0])*drdt(i+i0)*aa;
      // coefc[i]=id*h*(en-v[i+i0])*drdt(i+i0)*aa;
      // coefd[i]=-id*h*(ga-ka)*dror(i+i0);
      double dror = drdt[i+i0]/r[i+i0];
      coefa.push_back(-id*h*(ga+ka)*dror);
      coefb.push_back(-id*h*(en+2*c2-v[i+i0])*drdt[i+i0]*alpha);
      coefc.push_back(id*h*(en-v[i+i0])*drdt[i+i0]*alpha);
      coefd.push_back(-id*h*(ga-ka)*dror);
      for (int j=0; j<AMO; j++){
        em[i][j]=ie[i][j];
      }
      em[i][i]=em[i][i]-coefd[i];
    }


    // //inverts the matrix!  invfm = Inv(fm)
    // double invem[amo2][amo2]={0};
    // //invertmat(em,invem,amo2);
    // invertMatrix(em,invem);
    // // XXX update to use GSL ? XXX
    MAT_invertMatrix(em); //from here on, em is the inverted matrix!


    std::vector<double> s(AMO);
    std::vector< std::vector<double> > fm(AMO,std::vector<double>(AMO));
    for (int i=0; i<AMO; i++){
      s[i]=-ia[i]*u0;
      for (int j=0; j<AMO; j++){
        fm[i][j]=ie[i][j]-coefb[i]*em[i][j]*coefc[j];
        s[i]=s[i]-coefb[i]*em[i][j]*ia[j]*v0;
      }
      fm[i][i]=fm[i][i]-coefa[i];
    }


    //inverts the matrix!  invfm = Inv(fm)
    // XXX use vector!
    //double invfm[amo2][amo2]={0};
    //invertmat(fm,invfm,amo2);
    //invertMatrix(fm,invfm);
    // XXX update to use GSL ? XXX
    MAT_invertMatrix(fm); //from here on, fm is the inverted matrix!


    //writes u(r) in terms of coefs and the inverse of fm
    // P(r) = r^gamma u(r)
    // double us[amo2];
    std::vector<double> us(AMO);
    for (int i=0; i<AMO; i++){
      us[i]=0;
      for (int j=0; j<AMO; j++){
        us[i]=us[i]+fm[i][j]*s[j];
      }
    }


    //writes v(r) in terms of coefs + u(r)
    // Q(r) = r^gamma v(r)
    // double vs[amo2];
    std::vector<double> vs(AMO);
    for (int i=0; i<AMO; i++){
      vs[i]=0;
      for (int j=0; j<AMO; j++){
        vs[i]=vs[i]+em[i][j]*(coefc[j]*us[j]-ia[j]*v0); //XXX
        //XXX here: is this the large cancellation?
      }
    }


    //writes wavefunction: P= r^gamma u(r) etc..
    for (int i=0; i<AMO; i++){
      p[i+i0]=pow(r[i+i0],ga)*us[i];
      q[i+i0]=pow(r[i+i0],ga)*vs[i];
      //XXX r: to array!
    }

    //re-sets 'starting point' for next ln
    u0=us[AMO-1];
    v0=vs[AMO-1];

  }// END for (int ln=0; ln<nol; ln++)  [loop through outint `nol' times]


  // calls adamsmoulton to finish integration from (nol*AMO+1) to ctp
  int na=nol*AMO+1;
  if (ctp>na){
    adamsMoulton(p,q,en,v,ka,r,drdt,h,NGP,na,ctp,alpha);
  }


  return 0;
}  // END outint






//******************************************************************
int inwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int ctp, int pinf, double alpha)
/*
Program to start the INWARD integration.
Starts from Pinf, and uses an expansion(?) to go to (pinf-AMO)
Then, it then call ADAMS-MOULTON, to finished (from nol*AMO+1 to ctp)
*/
{

  const int nx=30;	//order of the expansion coeficients in 'inint'  (15 orig.)
  const float nxepsp=1e-18; // PRIMARY convergance for expansion in `inint'' (10^-8)
  const float nxepss=1e-3;  // SECONDARY convergance for expansion in `inint'' (10e-3)
  // XXX treat as variable here? Param in the .h ??

  double alpha2 = pow(alpha,2);
  double cc = 1./alpha;
  double c2 = 1./alpha2;

  double lambda=sqrt(-en*(2+en*alpha2)); //XXX alpha^2 ??
  double zeta=-v[pinf]*r[pinf];
  double sigma=(1+en*alpha2)*(zeta/lambda);
  double Ren=en+c2;   //total relativistic energy

  //Generates the expansion coeficients for asymptotic wf
  // up to order NX (nx is 'param')
  // double bx[nx];
  // double ax[nx];
  std::vector<double> bx(nx); //nx = ?? from 'params'?
  std::vector<double> ax(nx);
  bx[0]=(ka+(zeta/lambda))*(alpha/2);
  for (int i=0; i<nx; i++){
    ax[i]=(ka+(i+1-sigma)*Ren*alpha2-zeta*lambda*alpha2)*bx[i]*cc/((i+1)*lambda);
    if (i<(nx-1)){
      bx[i+1]=(pow(ka,2)-pow((i+1-sigma),2)-pow(zeta,2)*alpha2)
              *bx[i]/(2*(i+1)*lambda);
    }
  }


  //Generates last `AMO' points for P and Q [actually AMO+1?]
  double f1=sqrt(1.+en*alpha2/2.);
  double f2=sqrt(-en/2.)*alpha;
  for (int i=pinf; i>=(pinf-AMO); i--){      // double check end point!
    double rfac=pow(r[i],sigma)*exp(-lambda*r[i]);
    double ps=1.;
    double qs=0.;
    double rk=1.;
    for (int k=0; k<nx; k++){    //this will loop until a) converge, b) k=nx
      rk=rk*r[i]; //XXX replace w/ array!
      ps=ps+(ax[k]/rk);
      qs=qs+(bx[k]/rk);
      double xe=fmax(fabs((ax[k]/rk)/ps),fabs((bx[k]/rk)/qs));
      if (xe<nxepsp){
        k=nx;          //reached convergance
      }
      else if (k==(nx-1)){
        if (xe>nxepss){
          printf("WARNING: Asymp. expansion in ININT didn't converge:"
                 " %i, %i, %.2e\n"
            ,i,k,xe);
        }
      }
    }
    p[i]=rfac*(f1*ps+f2*qs);
    q[i]=rfac*(f2*ps-f1*qs);  //XXX here? the 'small cancellation'(?)
  }


  //calls adams-moulton
  if ((pinf-AMO-1)>=ctp){
    adamsMoulton(p,q,en,v,ka,r,drdt,h,NGP,pinf-AMO-1,ctp,alpha);
  }

  return 0;
}  // END inint


//******************************************************************************
int adamsMoulton(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int ni, int nf, double alpha)
/*
//program finishes the INWARD/OUTWARD integrations (ADAMS-MOULTON)
  //- ni is starting (initial) point for integration
  //- nf is end (final) point for integration (nf=ctp)
*/
{

  double c2 = 1./pow(alpha,2); // c^2 - just to shorten code

  // double ama[amo2];    //AM coefs.
  std::vector<double> ama(AMO);
  double amd,amaa;
  getAdamsCoefs(ama,amd,amaa);  // loads the coeficients! //XXX update!

  //this just checks that all working..prints out coefs.
  if (dodebug){
    for (int i=0; i<AMO; i++)
      printf("AMA[%i]=%f\n",i,ama[i]);
    printf("amd=%.0f, amaa=%.0f\n",amd,amaa);
  }

  int inc;      //'increment' for integration (+1 for forward, -1 for backward)
  int nosteps;    // number of steps integration should make
  if (nf>ni){
    inc=1;
    nosteps=nf-ni+1;     //check the "+1"...
  }
  else if (nf<ni){
    inc=-1;
    nosteps=ni-nf+1;
  }
  else{
    printf("WARNING: ni=nf in adamsmoulton.. no further integration");
    return 1;
  }


  // double dp[NGP],dq[NGP];      //create arrays for wf derivatives
  // double amcoef[amo2];
  //create arrays for wf derivatives
  std::vector<double> dp(NGP),dq(NGP); //XXX is this syntax valid?
  std::vector<double> amcoef(AMO);
  int k1=ni-inc*(AMO);
  for (int i=0; i<AMO; i++){//nb: k1 is iterated
    double dror = drdt[k1]/r[k1];
    dp[i]=inc*(-ka*dror*p[k1]-alpha*((en+2*c2)-v[k1])*drdt[k1]*q[k1]);
    dq[i]=inc*(ka*dror*q[k1]+alpha*(en-v[k1])*drdt[k1]*p[k1]);
    //XXX drdt etc. update!
    amcoef[i]=(h/amd)*ama[i];
    k1+=inc;
  }


  //integrates the function from ni to the c.t.p
  double a0=(h/amd)*amaa;
  int k2=ni;
  for (int i=0; i<nosteps; i++){    //double check! - end point should be ctp! (inclusive)
    double dror = drdt[k2]/r[k2];
    double dai=-inc*(ka*dror);
    double dbi=-inc*alpha*(en+2*c2-v[k2])*drdt[k2];
    double dci=inc*alpha*(en-v[k2])*drdt[k2];
    double ddi=-dai;
    double det = 1-a0*a0*(dbi*dci-dai*ddi);
    double sp=p[k2-inc];
    double sq=q[k2-inc];
    for (int l=0; l<AMO; l++){
      sp=sp+amcoef[l]*dp[l];
      sq=sq+amcoef[l]*dq[l];
    }
    p[k2]=(sp+a0*(dbi*sq-ddi*sp))/det;
    q[k2]=(sq+a0*(dci*sp-dai*sq))/det;
    for (int l=0; l<(AMO-1); l++){    //loads next 'first' k values (?)
      dp[l]=dp[l+1];
      dq[l]=dq[l+1];
    }
    dp[AMO-1]=dai*p[k2]+dbi*q[k2];    //loads next 'first' deriv's (?)
    dq[AMO-1]=dci*p[k2]+ddi*q[k2];
    k2+=inc;
  }

  return 0;
}  // END adamsmoulton



//**************************************************************************
//**************************************************************************
//**************************************************************************
//**************************************************************************



//       Numerical DE solve by integration (adams) coeficients


//******************************************************************
// coeficients for the ADAMS-MOULTON routine
int getAdamsCoefs(std::vector<double> &mia, double &mid, double &miaa)
{

// coefs for Adams..

//XXX XXX XXX Make array of MAX domension! then, just 'input' the correct #

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
}  // END AMcoefs



//******************************************************************
// coeficients for the OUTINT routine
//int OIcoefs(double (*oie)[amo2], double *oia, double &oid)
int getOutwardCoefs(std::vector< std::vector<double> > &oie,
    std::vector<double> &oia, double &oid)
{

  //XXX XXX XXX Make array of MAX domension! then, just 'input' the correct #
  // ?? little harder here..

// coefs for outint..
  //int AMO = oia.size();
  if (AMO==8){
    double tie[8][8]={-1338,2940,-2940,2450,-1470,588,-140,15,
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
      //oia.push_back(tia[i]);
      for (int j=0; j<AMO; j++){
        oie[i][j]=tie[i][j];
      }
    }
  }
  else if (AMO==7){
    double tie[7][7]={-609,1260,-1050,700,-315,84,-10,
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
    double tie[6][6]={-77,150,-100,50,-15,2,
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
    double tie[5][5]={-65,120,-60,20,-3,
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
}  //   END OIcoefs