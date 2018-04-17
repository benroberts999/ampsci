#include "HF_hartree.h"


//******************************************************************************
int HF_hartreeCore(ElectronOrbitals &wf, double eps_hartree)
{


  //Fill the electron part of the potential, using Greens PRM for initial approx
  double Gh,Gd;  //Green potential parameters
  PRM_defaultGreen(wf.Z,Gh,Gd); //Get default values for Green potential
  for(int i=0; i<wf.ngp; i++) wf.vdir.push_back(PRM_green(wf.Z,wf.r[i],Gh,Gd));

  //First step: Solve each core state using parameteric potential
  wf.solveInitialCore(1);

  //Hartree loop:
  int num_its=0;
  for(int n=0; n<max_hartree; n++){

    double eta=0.50;

    std::vector<double> vdir_old = wf.vdir;
    std::vector<double> vdir_new;
    formNewVdir(wf,vdir_new);
    for(int j=0; j<wf.ngp; j++){
      wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];
    }

    double prev_e = 0;
    for(size_t i=0; i<wf.nlist.size(); i++) prev_e += wf.en[i]/wf.nlist.size();

    for(size_t i=0; i<wf.nlist.size(); i++){
      double del_e=0;
      for(int j=0; j<wf.ngp; j++)
        del_e += (wf.vdir[j]-vdir_old[j])*
        (pow(wf.p[i][j],2) + pow(wf.q[i][j],2))*wf.drdt[j];
      del_e*=wf.h;
      double new_e = wf.en[i] + 1*del_e;
      if(new_e>0)new_e=-0.1;
      wf.reSolveLocalDirac(i,new_e,3); //only go to 1/10^3 - do better at end!
    }

    double next_e = 0;
    for(size_t i=0; i<wf.nlist.size(); i++) next_e += wf.en[i]/wf.nlist.size();

    double delta_hartree = (next_e-prev_e)/(next_e*eta);
    num_its = n+1;
    printf("Hart it:%3i,  del=%6.0e\n",num_its,delta_hartree);

    if(fabs(delta_hartree)<eps_hartree) break;
  }

  //re-run solve Dirac to higher convergance level after Hart pot. ok
  for(size_t i=0; i<wf.nlist.size(); i++) wf.reSolveLocalDirac(i,0,14);

  formNewVdir(wf,wf.vdir,false);


  std::ofstream ofile;
  ofile.open("pot.txt");
  ofile<<"r Gr Vh Z/r 1/r\n";
  for(int i=0; i<wf.ngp; i++){
    ofile<<wf.r[i]<<" "
      <<-PRM_green(wf.Z,wf.r[i],Gh,Gd)-wf.vnuc[i]<<" "
      <<-wf.vdir[i]-wf.vnuc[i]<<" "
      <<wf.Z/wf.r[i]<<" "<<1./wf.r[i]<<"\n";

  }
  ofile.close();


  return num_its;
}



//******************************************************************************
int formNewVdir(ElectronOrbitals wf, std::vector<double> &vdir_new, bool core)
//XXX OK? lots memory? Reference?
//XXX core list!!! make part of class!
{

  vdir_new.clear();
  vdir_new.resize(wf.ngp);

  //XXX use core_list instead!!! XXX ??
  int Ncore=0;
  for(size_t i=0; i<wf.nlist.size(); i++){
    int ka = wf.klist[i];
    int twoj = 2*abs(ka)-1;
    Ncore += twoj+1;
  }


  // //a=-1 means assume vdir same for all orbitals!
  double f=1;
  if(core) f = 1. - (1.)/Ncore;

  std::vector<double> rho(wf.ngp);
  for(size_t i=0; i<wf.nlist.size(); i++){
    int ka = wf.klist[i];
    int twoj = 2*abs(ka)-1;
    for(int j=0; j<wf.ngp; j++){
      rho[j] += (twoj+1)*(pow(wf.p[i][j],2) + pow(wf.q[i][j],2));
      //XXX assumes closed shell!? ia!
    }
  }

  for(int j=0; j<wf.ngp; j++){
    double r = wf.r[j];
    double v_tmp = 0;
    for(int k=0; k<wf.ngp; k++){
      double rp = wf.r[k];
      double rm = std::max(r,rp);//can be little more clever, slight speedup
      v_tmp += (rho[k]/rm)*wf.drdt[k];
    }
    vdir_new[j] = f*v_tmp*wf.h;
  }

  return 0;
}
