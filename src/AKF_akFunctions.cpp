#include "AKF_akFunctions.h"
namespace AKF{

//******************************************************************************
double CLkk(int L, int ka, int kb)
{
  int la = ATI::l_k(ka);
  int lb = ATI::l_k(kb);
  int two_ja = ATI::twoj_k(ka);
  int two_jb = ATI::twoj_k(kb);
  double ja = 0.5*two_ja;
  double jb = 0.5*two_jb;

  if((la+lb+L)%2!=0) return 0; //Parity rule
  if((la+lb<L)||(abs(la-lb)>L)) return 0; //triangle rule (l)

  double tjB = WIG::threej(jb,ja,L,-0.5,0.5,0);
  return (2*ja+1)*(2*jb+1)*(2*L+1)*pow(tjB,2);
}


//******************************************************************************
void writeToTextFile(
  std::string fname,
  std::vector< std::vector< std::vector<float> > > &AK,
  std::vector<std::string> nklst,
  double qmin, double qmax,
  double demin, double demax
)
{
  int desteps = AK.size();   //dE
  int num_states = AK[0].size();  //nk
  int qsteps = AK[0][0].size();//q

  double qMeV = (1.e6/(FPC::Hartree_eV*FPC::c));
  double keV = (1.e3/FPC::Hartree_eV);

  std::ofstream ofile;
  ofile.open(fname+".txt");
  ofile<<"dE(keV) q(MeV) ";
  for(size_t i=0; i<nklst.size(); i++) ofile<<nklst[i]<<" ";
  ofile<<"\n\n";
  for(int i=0; i<desteps; i++){
    for(int k=0; k<qsteps; k++){
      double x = double(k)/(qsteps-1);
      if(qsteps==1) x=0;
      double q = qmin*pow(qmax/qmin,x);
      double y = double(i)/(desteps-1);
      if(desteps==1) y=0;
      double dE = demin*pow(demax/demin,y);
      ofile<<dE/keV<<" "<<q/qMeV<<" ";
      for(int j=0; j<num_states; j++){
        ofile<<AK[i][j][k]<<" ";
      }
      ofile<<"\n";
    }
    if(qsteps>1)ofile<<"\n";
  }
  ofile.close();
}


//******************************************************************************
/*
Helper functions for the binary read in/out
*/
template<typename T>
int binary_rw(std::fstream& stream, T& value, bool write){
  if(write) stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
  else stream.read(reinterpret_cast<char*>(&value), sizeof(T));
  return 0;
}
template<typename T>
int binary_str_rw(std::fstream& stream, T& value, bool write){
  if(write){
    size_t temp_len = value.length();
    stream.write(reinterpret_cast<const char*>(&temp_len), sizeof(size_t));
    stream.write(value.c_str(), value.length());
  }else{
    size_t temp_len;
    stream.read(reinterpret_cast<char*>(&temp_len), sizeof(size_t));
    char* tvalue = new char[temp_len+1];
    stream.read(tvalue, temp_len);
    tvalue[temp_len] = '\0'; //null 'end of string' character
    value = tvalue;
    delete [] tvalue;
    return 0;
  }
  return 0;
}


//******************************************************************************
int akReadWrite(std::string fname, bool write,
  std::vector< std::vector< std::vector<float> > > &AK,
  std::vector<std::string> &nklst,
  double &qmin, double &qmax,
  double &dEmin, double &dEmax)
{
  std::fstream iof;
  fname=fname+".bin";
  if(write) iof.open(fname,std::ios_base::out|std::ios_base::binary);
  else      iof.open(fname,std::ios_base::in |std::ios_base::binary);

  if(write){
    int nde = AK.size();   //dE
    int ns = AK[0].size();  //nk
    int nq = AK[0][0].size();//q
    binary_rw(iof,nde,write);
    binary_rw(iof,ns,write);
    binary_rw(iof,nq,write);
  }else{
    int nq,ns,nde;
    binary_rw(iof,nde,write);
    binary_rw(iof,ns,write);
    binary_rw(iof,nq,write);
    AK.resize(nde,std::vector< std::vector<float> >(ns,std::vector<float>(nq)));
    nklst.resize(ns);
  }
  binary_rw(iof,qmin,write);
  binary_rw(iof,qmax,write);
  binary_rw(iof,dEmin,write);
  binary_rw(iof,dEmax,write);
  for(size_t ie=0; ie<AK.size(); ie++){
    for(size_t in=0; in<AK[0].size(); in++){
      if(ie==0) binary_str_rw(iof,nklst[in],write);
      for(size_t iq=0; iq<AK[0][0].size(); iq++){
        binary_rw(iof,AK[ie][in][iq],write);
      }
    }
  }

  return 0;
}


//******************************************************************************
int calculateK_nk(ElectronOrbitals &wf, int is, int max_L, double dE,
  std::vector< std::vector<std::vector<float> > > &jLqr_f,
  std::vector< std::vector<float> > &K_nk, double Zeff)
/*
Zeff is '-1' by default. If Zeff > 0, will solve w/ Zeff model
*/
{
  if (is>=(int)wf.p.size()) return 1; //should never occur, just safety?

  ContinuumOrbitals cntm(wf);// create cntm object [survives locally only]

  int k = wf.klist[is];
  int l = ATI::l_k(k);

  int qsteps = (int)jLqr_f[0].size();

  //Calculate continuum wavefunctions
  double ec = dE+wf.en[is];
  cntm.clear();
  int lc_max = l + max_L;
  int lc_min = l - max_L;
  if(lc_min<0) lc_min = 0;
  if(ec>0){
    if(Zeff>0) cntm.solveZeffContinuum(ec,Zeff,lc_min,lc_max); //Zeff version
    else       cntm.solveLocalContinuum(ec,lc_min,lc_max);
  }

  // Generate AK for each L, lc, and q
  // L and lc are summed, not stored indevidually
  std::vector<float> AK_nk_q(qsteps);
  for(int L=0; L<=max_L; L++){
    for(size_t ic=0; ic<cntm.klist.size(); ic++){
      int kc = cntm.klist[ic];
      double dC_Lkk = CLkk(L,k,kc); //XXX new formula!
      if(dC_Lkk==0) continue;
      #pragma omp parallel for
      for(int iq=0; iq<qsteps; iq++){
        double a = 0;
        double jLqr = 0;
        if(cntm.p.size()>0){
          if(ec<=0) std::cout<<"ERROR 244: !?!?\n";
          int maxj = wf.pinflist[is]; //don't bother going further
          //Do the radial integral:
          a=0;
          for(int j=0; j<maxj; j++){
            jLqr = jLqr_f[L][iq][j];
            a += (wf.p[is][j]*cntm.p[ic][j] + wf.q[is][j]*cntm.q[ic][j])
                 *jLqr*wf.drdt[j];// *h below!
          }
        }
        AK_nk_q[iq] += dC_Lkk*pow(a*wf.h,2);
      } //q
    } // END loop over cntm states (ic)
  } // end L loop
  K_nk.push_back(AK_nk_q);
  cntm.clear(); //deletes cntm wfs for this energy
  return 0;
}


//******************************************************************************
int calculateKpw_nk(ElectronOrbitals &wf, int nk, double dE,
  std::vector< std::vector<float> > &jl_qr,
  std::vector< std::vector<float> > &K_nk
)
/*
For plane-wave final state.
Only has f-part....Can restore g-part, but need to be sure of plane-wave!
Chi(q) - Int[ f_nk*j_l(qr)*r , {r,0,inf}]
Should be called once per initial state

XXX Note sure if correct! esp, (q) angular part!? XXX

*/
{
  if (nk>=(int)wf.p.size()) return 1; //should never occur

  int kappa = wf.klist[nk];
  int twoj = ATI::twoj_k(kappa);

  int qsteps = (int)jl_qr.size();
  std::vector<float> tmpK_q(qsteps);

  double eps = dE - wf.en[nk];
  int maxir = wf.pinflist[nk]; //don't bother going further

  for(int iq=0; iq<qsteps; iq++){
    if(eps<=0) break;
    double chi_q=0;
    for(int ir=0; ir<maxir; ir++){
      chi_q += wf.p[nk][ir]*jl_qr[iq][ir]*wf.r[ir]*wf.drdt[ir];
    }
    chi_q *= wf.h;
    tmpK_q[iq] = (2./M_PI)*(twoj+1)*pow(chi_q,2)*sqrt(2.*eps);
    //tmpK_q[iq] = pow(4*3.14159,2)*pow(chi_q,2); //XXX XXX just cf KOPP
  }

  K_nk.push_back(tmpK_q);
  return 0;
}


//******************************************************************************
void sphericalBesselTable(
  std::vector< std::vector< std::vector<float> > > &jLqr_f,
  int max_L,
  double qmin, double qmax, int qsteps,
  std::vector<double> &r)
{
  std::cout<<std::endl;
  int ngp = (int)r.size();
  jLqr_f.resize(max_L+1, std::vector< std::vector<float> >
    (qsteps, std::vector<float>(ngp)));
  for(int L=0; L<=max_L; L++){
    std::cout<<"\rCalculating spherical Bessel look-up table for L="
    <<L<<"/"<<max_L<<" .. "<<std::flush;
    #pragma omp parallel for
    for(int iq=0; iq<qsteps; iq++){
      double x=double(iq)/(qsteps-1);
      double q = qmin*pow(qmax/qmin,x);
      for(int ir=0; ir<ngp; ir++){
        double tmp = gsl_sf_bessel_jl(L, q*r[ir]);
        // If q(dr) is too large, "missing" j_L oscillations
        //(overstepping them). This helps to fix that.
        //By averaging the J_L function. Note: only works if wf is smooth
        int num_extra = 0;
        if(ir<ngp-1){
          double qdrop = q*(r[ir+1]-r[ir])/M_PI;
          double min_qdrop = 0.05; // require 20 pts per half wavelength!
          if(qdrop>min_qdrop) num_extra = int(qdrop/min_qdrop)+1;
        }
        {//Include 'extra' points into j_L (avg):
          for(int i=0; i<num_extra; i++){
            double b = (i+1.)/(num_extra+1.);
            double a=1.-b;
            double qrtmp = q*(a*r[ir]+b*r[ir+1]);
            tmp += gsl_sf_bessel_jl(L, qrtmp);
          }
          tmp /= (num_extra+1);
        }

        jLqr_f[L][iq][ir] = tmp;
      }
    }
  }
  std::cout<<"done\n";
}


}//End namespace
