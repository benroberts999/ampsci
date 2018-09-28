#include "ADAMS_solveLocalBS.h"
/*
Thu 06 Nov 2014 23:29:57 AEDT. Updated 2018 to go past ctp.
Program to solve single-electron bound-state Dirac problem for a (given)
local, central potential.
Based on method presented in book by W. Johnson.
Employs the Adams-Moulton method.
solveDBS is the main routine that is called from elsewhere.
All other functions called by solveDBS.

===== To do :: Jan 2018 =====
  * p,q -> f,g! Check the q cancellation!

*/

namespace ADAMS{

bool debug=false; //if true, will print progress messages.
//parameter: Adams-Moulton ''order'' (number of points)



//******************************************************************************
int solveDBS(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, double Z, int n, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int &pinf, int &its, double &eps, double alpha, int log_dele_or)
/*
Solves local, spherical bound state dirac equation using Adams-Moulton method.
Based on method presented in book by W. R. Johnson:
  W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007)
I have added a few extensions to this method. In particular, I integrate past
the classical turning point. (See below)

See also:
 * https://en.wikipedia.org/wiki/Linear_multistep_method
 * Hairer, Ernst; Nørsett, Syvert Paul; Wanner, Gerhard (1993),
   Solving ordinary differential equations I: Nonstiff problems (2nd ed.),
   Berlin: Springer Verlag, ISBN 978-3-540-56670-0.
 * Quarteroni, Alfio; Sacco, Riccardo; Saleri, Fausto (2000),
   Matematica Numerica, Springer Verlag, ISBN 978-88-470-0077-3.
 * http://mathworld.wolfram.com/AdamsMethod.html

Rough description of method:
1. Start with initial 'guess' of energy
2. Find the "practical infinity" (psi~0), and the Classical turning point [e=v]
3. Performs 'inward' integration (Adams Moulton). Integrates from the practical
   infinity inwards to d_ctp points past the classical turning point (ctp).
4. Performs 'outward' integration (Adams Moulton). Integrates from 0
   outwards to d_ctp points past the ctp.
5. Matches the two functions around ctp, by re-scaling the 'inward' solution
   (for P). Uses a weighted mean, with weights given by distance from ctp.
6. Checks the number of nodes the wf has. If too many or too few nodes, makes
   a large change to the energy and tries again (from step 2).
   If the correct number of nodes, uses perturbation theory to make minor
   corrections to the energy to 'zoom in' (matching the in/out solution for q),
   then re-starts from step 2.
Continues until this energy adjustment falls below a prescribed threshold.

There are two energy convergence parameters, primary and secondary.
If cannot converge to primary value, but better than secondary, will not give
an error.

*/
{
  //Parameters. Put somwhere else?
  double delep=1e-15; //PRIMARY convergence parameter [energy] (10^-11)
  double deles=1e-6; //SECONDAY convergence parameter [energy]	(X)
  const int ntry=32;        // Max # attempts at converging [sove bs] (30)
  const double alr=800;     // ''assymptotically large r [kinda..]''  (=800)
  const double lde=0.3;     // 'large' energy variations (0.1 => 10%)
  const int d_ctp=5;        //Num points past ctp +/- d_ctp.

  if(log_dele_or>0){
    delep=1./pow(10,log_dele_or);
    deles=10*delep;
  }

  //Checks to see if legal n is requested. If not, increases n, re-calls
  //Should I do this? If there's a logic problem, probably better to know?
  if ( (fabs(ka)<=n) && (ka!=n) ){
    if(debug) printf("\nRunning SolveDBS for state %i %i:\n",n,ka);
  }else{
    printf("\nSate %i %i does not exist..\n",n,ka);
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
  int ctp=-1;                 //classical turning point
  bool converged=false;
  double delta_en=0;
  its=0;            //numer of iterations (for this n,ka)

  while (!converged){

    //Find the practical infinity 'pinf'
    //Step backwards from the last point (NGP-1) until
    // (V(r) - E)*r^2 >  alr    (alr = "asymptotically large r")
    pinf=NGP-1;
    while((en-v[pinf])*pow(r[pinf],2) + alr < 0) pinf--; //pinf=pinf-1;
    if((pinf==NGP-1)&& debug)
      printf("WARNING: pract. inf. = size of box for %i %i\n",n,ka);
    if(debug)
      printf("Practical infinity (i=%i): Pinf=%.1f a.u.\n",pinf,r[pinf]);

    //Find classical turning point 'ctp'
    //Step backwards from the "practical infinity" until
    //  V(r) > E        [nb: both V and E are <0]
    ctp=pinf;
    while ( (en-v[ctp]) < 0 ){
      ctp--; //=ctp-1;
      if (ctp-d_ctp<=0){
        //fails if ctp<0, (or ctp>pinf?)
        printf("FAILURE 96 in solveDBS: No classical region?\n");
        return 1;
      }
    }
    if(ctp+d_ctp>=pinf){
      //Didn't find ctp! Does this ever happen? Yes, if energy guess too wrong
      if(debug)printf("FAILURE: Turning point at or after pract. inf. \n");
      if(debug)printf("ctp=%i, d_ctp=%i, pinf=%i, NGP=%i\n",ctp,d_ctp,pinf,NGP);
      return 1;
    }
    if(debug) printf("Classical turning point (i=%i): ctp=%.1f a.u.\n",
                ctp,r[ctp]);
    if(debug) printf("%i %i: Pinf= %.1f,  en= %f\n",n,ka,r[pinf],en);

    //Temporary vectors for in/out integrations
    std::vector<double> pin(NGP),qin(NGP),pout(NGP),qout(NGP);
    //Perform the "inwards integration":
    inwardAM(pin,qin,en,v,ka,r,drdt,h,NGP,ctp-d_ctp,pinf,alpha);
    //Perform the "outwards integration"
    outwardAM(pout,qout,en,v,Z,ka,r,drdt,h,NGP,ctp+d_ctp,alpha);

    //Find the re-scaling factor using a weighted average:
    double rescale = pout[ctp]/pin[ctp];
    double denom = 1;
    for(int i=1; i<=d_ctp; i++){
      //re-scale from weigted average. Weight = distance from ctp
      rescale += 0.5*(pout[ctp+i]/pin[ctp+i] + pout[ctp-i]/pin[ctp-i])/(i+1);
      denom += 1./(i+1);
    }
    rescale /= denom;
    if(debug) std::cout<<"denom="<<denom<<" rescale="<<rescale<<"\n";

    //re-scale the inward solution to match the outward solution (for P):
    for(int i=ctp-d_ctp; i<=pinf; i++){
      pin[i]*=rescale;
      qin[i]*=rescale;
    }

    //Join the in and outward solutions. "Meshed" around ctp +/- d_ctp
    for(int i=0; i<ctp-d_ctp; i++){
      p[i] = pout[i];
      q[i] = qout[i];
    }
    for(int i=ctp+d_ctp+1; i<=pinf; i++){
      p[i] = pin[i];
      q[i] = qin[i];
    }
    for(int i=ctp-d_ctp; i<=ctp+d_ctp; i++){
      //"Mesh" in the intermediate region
      double B;
      if(d_ctp==0) B=0.5;
      else B = (double((i-ctp)+d_ctp))/(2.*d_ctp);
      double A = 1-B;
      p[i] = A*pout[i] + B*pin[i];
      q[i] = A*qout[i] + B*qin[i];
    }

    //Count the number of nodes (zeros) the wf has.
    //Just counts the number of times wf changes sign
    int nozeros=0;
    double sp=p[1];
    double spn;
    for(int i=2; i<pinf; i++){
      spn=p[i];
      if (sp*spn<0) nozeros++;
      if(spn!=0)    sp=spn;
    }
    if(debug) printf("Nodes=%i\n",nozeros);

    //checks to see if there are too many/too few nodes.
    //Makes large adjustment to energy until in correct region
    //Then, once in 'correct' region, uses perturbation theory to make small
    // changes to the energy
    double etemp;
    if(nozeros>inodes){
      //Too many nodes:
      more++;
      if((more==1)||(en<eupper)) eupper=en;
      etemp=(1+lde)*en;
      if((less!=0)&&(etemp<elower)) etemp=0.5*(eupper+elower);
      delta_en=fabs((en-etemp)/en);
      en=etemp;
    }else if (nozeros<inodes){
      //too few nodes:
      less++;
      if((less==1)||(en>elower)) elower=en;
      etemp=(1-lde)*en;
      if((more!=0)&&(etemp>eupper)) etemp=0.5*(eupper+elower);
      delta_en=fabs((en-etemp)/en);
      en=etemp;
    }else{
      // XXX Maybe, put this into a function!
      // correct number of nodes.
      //From here, use perturbation theory to fine-time the energy
      if(debug) printf("Correct number of nodes, starting P.T.\n");
      std::vector<double> ppqq(NGP);
      for (int i=0; i<=pinf; i++){
        ppqq[i]=p[i]*p[i]+q[i]*q[i];    // XXX add alpha here if need!
      }
      anorm=INT::integrate(ppqq,drdt,h,0,pinf);
      if(debug) printf("anrom=%.5f\n",anorm);
      //Use perturbation theory to work out delta En
      // delta E = c*P(r)*[Qin(r)-Qout(r)]
      // nb: wf not yet normalised!
      double p_del_q = p[ctp]*(qin[ctp]-qout[ctp]);
      double denom=1;
      //weighted average around ctp:
      for(int i=1; i<=d_ctp; i++){
        p_del_q += 0.5*(
                    p[ctp+i]*(qin[ctp+i]-qout[ctp+i]) +
                    p[ctp-i]*(qin[ctp-i]-qout[ctp-i])
                  )/(i+1.);
        denom += 1./(i+1.);
      }
      p_del_q /= denom;
      double de=  (1./alpha) * p_del_q / anorm ;
      delta_en=fabs(de/en);
      etemp = en + de;
      if(debug) printf("de=%.3e, en=%.5f, et=%.5f, el=%.5f, e=%.5f\n",
                de,en,etemp,elower,eupper);
      if((less!=0)&&(etemp<elower)){
        delta_en=fabs((en-0.5*(en+elower))/en);
        en=0.5*(en+elower);
      }else if((more!=0)&&(etemp>eupper)){
        delta_en=fabs((en-0.5*(en+eupper))/en); //should be fabs? wasn't??
        en=0.5*(en+eupper);
      }else if(delta_en<delep){
        en=etemp;
        eps=delta_en;
        converged=true;
      }else{
        eps=delta_en;
        en=etemp;
      }
    }// END: if(nozeros>inodes

    its++; //increment 'number of iterations' counter
    if(debug) printf("Itteration number %i,  en= %f\n",its,en);

    if(its>=ntry){
      if(delta_en<deles){
        if(debug) printf("Wavefunction %i %i didn't fully converge after "
                    "%i iterations, but OK.\n",n,ka,ntry);
        eps=delta_en;
        converged=true; //Converged to "secondary" level, deles
      }else{
        printf("Wavefunction %i %i didn't converge after %i itterations.\n",n,
               ka,ntry);
        if(debug) printf("%i %i: Pinf= %.1f,  en= %f\n",n,ka,r[pinf],en);
        eps=delta_en;
        return 2;
      }
    }

  }// END: while (not converged)

  //normalises the wavefunction
  double an= 1/sqrt(anorm);
  if(debug) printf("An=%.5e\n",an);
  for(int i=0; i<=pinf; i++){
    p[i]=an*p[i];
    q[i]=an*q[i];
  }
  for(int i=pinf+1; i<NGP; i++){
    // kills remainders (just for safety)
    p[i]=0;
    q[i]=0;
  }

  if(debug){
    printf("NGP=%i, Size of box: Rmax=%.1f a.u., h=%f\n",NGP,r[NGP-1],h);
    printf("Converged for %i %i to %.3e after %i iterations;\n",n,ka,eps,its);
    printf("%i %i: Pinf(%i)=%.1f,  \nMy energy:  en= %.15f\n",n,ka,pinf,
           r[pinf],en);
  }

  return 0; // XXX XXX XXX have a code!!!!!
}

//*********************************************************
//*********************************************************
//*********************************************************

//******************************************************************************
int outwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, double Z, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int nf, double alpha)
/*
Program to start the OUTWARD integration.
Starts from 0, and uses an expansion(?) to go to (nol*AMO).
Then, it then call ADAMS-MOULTON, to finish (from nol*AMO+1 to nf = ctp+d_ctp)
*/
{
  const int nol=1; // # of outdir runs [finds first nol*AMO+1 points (3)]

  double az = Z*alpha;
  double c2 = 1./pow(alpha,2);
  double ga = sqrt(pow(ka,2)-pow(az,2));

  //initial wf values
  // P(r) = r^gamma u(r)
  // Q(r) = r^gamma v(r)
  double u0=1;
  double v0;
  if (ka>0) v0 = -(ga+ka)/az;
  else      v0 = az/(ga-ka);
  p[0]=0;
  q[0]=0;

  //Coeficients used by the method
  std::vector< std::vector<double> > ie(AMO,std::vector<double>(AMO));
  std::vector<double> ia(AMO);
  double id;
  getOutwardCoefs(ie,ia,id);

  // loop through and find first nol*AMO points of wf
  for (int ln=0; ln<nol; ln++){
    int i0=ln*AMO+1;

    //defines/populates em coefs
    std::vector<double> coefa,coefb,coefc,coefd;
    std::vector< std::vector<double> > em(AMO,std::vector<double>(AMO));
    for (int i=0; i<AMO; i++){
      double dror = drdt[i+i0]/r[i+i0];
      coefa.push_back(-id*h*(ga+ka)*dror);
      coefb.push_back(-id*h*(en+2*c2-v[i+i0])*drdt[i+i0]*alpha);
      coefc.push_back(id*h*(en-v[i+i0])*drdt[i+i0]*alpha);
      coefd.push_back(-id*h*(ga-ka)*dror);
      for (int j=0; j<AMO; j++) em[i][j]=ie[i][j];
      em[i][i]=em[i][i]-coefd[i];
    }
    // //inverts the em matrix
    MAT::invertMatrix(em); //from here on, em is the inverted matrix

    //defines/populates fm, s coefs
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
    //inverts the matrix!  fm =-> Inv(fm)
    MAT::invertMatrix(fm); //from here on, fm is the inverted matrix

    //writes u(r) in terms of coefs and the inverse of fm
    // P(r) = r^gamma u(r)
    std::vector<double> us(AMO);
    for (int i=0; i<AMO; i++){
      us[i]=0;
      for (int j=0; j<AMO; j++){
        us[i]=us[i]+fm[i][j]*s[j];
      }
    }

    //writes v(r) in terms of coefs + u(r)
    // Q(r) = r^gamma v(r)
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
    }

    //re-sets 'starting point' for next ln
    u0=us[AMO-1];
    v0=vs[AMO-1];

  }// END for (int ln=0; ln<nol; ln++)  [loop through outint `nol' times]

  //Call adamsmoulton to finish integration from (nol*AMO+1) to nf = ctp+d_ctp
  int na = nol*AMO+1;
  if (nf>na) adamsMoulton(p,q,en,v,ka,r,drdt,h,NGP,na,nf,alpha);

  return 0;
}


//******************************************************************
int inwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int nf, int pinf, double alpha)
/*
Program to start the INWARD integration.
Starts from Pinf, and uses an expansion(?) to go to (pinf-AMO)
Then, it then call ADAMS-MOULTON, to finish (from nol*AMO+1 to nf = ctp-d_ctp)
*/
{
  //order of the expansion coeficients in 'inint'  (15 orig.)
  const int nx=30;
  // PRIMARY convergance for expansion in `inint'' (10^-8)
  const float nxepsp=1e-18;
  // SECONDARY convergance for expansion in `inint'' (10e-3):
  const float nxepss=1e-3;

  double alpha2 = pow(alpha,2);
  double cc = 1./alpha;
  double c2 = 1./alpha2;

  double lambda=sqrt(-en*(2+en*alpha2));
  double zeta=-v[pinf]*r[pinf];
  double sigma=(1+en*alpha2)*(zeta/lambda);
  double Ren=en+c2;   //total relativistic energy

  //Generates the expansion coeficients for asymptotic wf
  // up to order NX (nx is 'param')
  std::vector<double> bx(nx);
  std::vector<double> ax(nx);
  bx[0]=(ka+(zeta/lambda))*(alpha/2);
  for(int i=0; i<nx; i++){
    ax[i]=(ka+(i+1-sigma)*Ren*alpha2-zeta*lambda*alpha2)*bx[i]*cc
         /((i+1)*lambda);
    if(i<(nx-1)) bx[i+1]=(pow(ka,2)-pow((i+1-sigma),2)-pow(zeta,2)*alpha2)
                *bx[i]/(2*(i+1)*lambda);
  }

  //Generates last `AMO' points for P and Q [actually AMO+1?]
  double f1=sqrt(1.+en*alpha2/2.);
  double f2=sqrt(-en/2.)*alpha;
  for (int i=pinf; i>=(pinf-AMO); i--){
    double rfac=pow(r[i],sigma)*exp(-lambda*r[i]);
    double ps=1.;
    double qs=0.;
    double rk=1.;
    for (int k=0; k<nx; k++){ //this will loop until a) converge, b) k=nx
      rk=rk*r[i];
      ps=ps+(ax[k]/rk);
      qs=qs+(bx[k]/rk);
      double xe=fmax(fabs((ax[k]/rk)/ps),fabs((bx[k]/rk)/qs));
      if (xe<nxepsp){
        k=nx; //reached convergance
      }
      else if (k==(nx-1)){
        if (xe>nxepss && debug)
          printf("WARNING: Asymp. expansion in ININT didn't converge:"
                 " %i, %i, %.2e\n",i,k,xe);
      }
    }
    p[i]=rfac*(f1*ps+f2*qs);
    q[i]=rfac*(f2*ps-f1*qs);  //XXX here? the 'small cancellation'(?)
  }

  //calls adams-moulton
  if ((pinf-AMO-1)>=nf){
    adamsMoulton(p,q,en,v,ka,r,drdt,h,NGP,pinf-AMO-1,nf,alpha);
  }

  return 0;
}


//******************************************************************************
int adamsMoulton(std::vector<double> &p, std::vector<double> &q, double &en,
    std::vector<double> v, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    int ni, int nf, double alpha)
/*
program finishes the INWARD/OUTWARD integrations (ADAMS-MOULTON)
  //- ni is starting (initial) point for integration
  //- nf is end (final) point for integration (nf=ctp+/-d_ctp)
*/
{
  double c2 = 1./pow(alpha,2); // c^2 - just to shorten code
  //AM coefs.
  std::vector<double> ama(AMO);
  double amd,amaa;
  getAdamsCoefs(ama,amd,amaa);  // loads the coeficients

  //this just checks that all working..prints out coefs.
  if (debug){
    for (int i=0; i<AMO; i++) printf("AMA[%i]=%f\n",i,ama[i]);
    printf("amd=%.0f, amaa=%.0f\n",amd,amaa);
  }

  int inc;     //'increment' for integration (+1 for forward, -1 for backward)
  int nosteps; // number of steps integration should make
  if (nf>ni){
    inc=1;
    nosteps=nf-ni+1;     //check the "+1"...
  }
  else if (nf<ni){
    inc=-1;
    nosteps=ni-nf+1;
  }
  else{
    printf("WARNING 611: ni=nf in adamsmoulton.. no further integration");
    return 1;
  }

  //create arrays for wf derivatives
  std::vector<double> dp(NGP),dq(NGP); //XXX is this syntax valid?
  std::vector<double> amcoef(AMO);
  int k1=ni-inc*(AMO);
  for (int i=0; i<AMO; i++){//nb: k1 is iterated
    double dror = drdt[k1]/r[k1];
    dp[i]=inc*(-ka*dror*p[k1]-alpha*((en+2*c2)-v[k1])*drdt[k1]*q[k1]);
    dq[i]=inc*(ka*dror*q[k1]+alpha*(en-v[k1])*drdt[k1]*p[k1]);
    amcoef[i]=(h/amd)*ama[i];
    k1+=inc;
  }

  //integrates the function from ni to the c.t.p
  double a0=(h/amd)*amaa;
  int k2=ni;
  for (int i=0; i<nosteps; i++){
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

//******************************************************************
int getAdamsCoefs(std::vector<double> &mia, double &mid, double &miaa)
/*
coeficients for the ADAMS-MOULTON routine
*/
{
//XXX XXX XXX Make array of MAX domension! then, just 'input' the correct #

  if (AMO==8){
    double tia[8]={-33953,312874,-1291214,3146338,-5033120,5595358,-4604594,
      4467094};
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
int getOutwardCoefs(std::vector< std::vector<double> > &oie,
    std::vector<double> &oia, double &oid)
/*
coeficients for the OUTINT routine
*/
{
//XXX XXX XXX Make array of MAX domension! then, just 'input' the correct #
// ?? little harder here..
//XX ? only down to 5?

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



}//end namespace
