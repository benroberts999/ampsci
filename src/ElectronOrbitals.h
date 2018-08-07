#ifndef _ORBITALS_H
#define _ORBITALS_H
#include <string>
#include <vector>
#include <cmath>
#include "FPC_physicalConstants.h"
#include "ATI_atomInfo.h"
#include "ADAMS_solveLocalBS.h"
// #include "adamsSolveLocalContinuum.h" //XXX sub-class???
#include <gsl/gsl_sf_fermi_dirac.h>
#include <algorithm> //for sort

const int NGP_DEFAULT=1000; //???

class ElectronOrbitals{

  public:

    ElectronOrbitals(int in_z, int in_a=0, int in_ngp=NGP_DEFAULT,
        double rmin=1e-6, double rmax=250., double var_alpha=1);
    ElectronOrbitals(std::string s_in_z, int in_a=0, int in_ngp=NGP_DEFAULT,
        double rmin=1e-6, double rmax=250., double var_alpha=1);

    std::vector< std::vector<double> > p;
    std::vector< std::vector<double> > q;
    std::vector<double> en;

    std::vector<double> r;
    std::vector<double> drdt;
    std::vector<double> dror;
    int ngp;
    double h;

    //std::vector<double> v;    //total ??
    std::vector<double> vnuc;
    std::vector<double> vdir; //direct/local part of the electron potential

    int Z,A;
    std::string atom;

    //double var_alpha; // like this?
    double alpha;

    int max_n;
    int max_l;
    int num_states; //?

    std::vector<int> core_list;
    std::vector<double> core_ocf; //occupancy fraction
    int num_core;

    std::vector<int> nlist;
    std::vector<int> klist;

    std::vector<int> pinflist;
    std::vector<int> itslist;
    std::vector<double> epslist;

    int hydrogenLike(int in_max_n, int in_max_l=100);

    int solveLocalDirac(int n, int k, double en_a, int log_dele_or=0);
    int reSolveLocalDirac(int i, double e_a=0, int log_dele_or=0);

    double diracen(double z, double n, int k);
    int sphericalNucleus(double rnuc=0);
    int fermiNucleus(double t=0, double c=0);

    int getRadialIndex(double r_target);

    int determineCore(std::vector<std::string> str_core);
    int solveInitialCore(int log_dele_or=0);
    bool isInCore(int n, int k);

    int sortedEnergyList(std::vector<int> &sort_list);


  private:

    //Number of grid points:


    int JohnsonRadialGrid(double r0=1.e-6, double rmax=250.);
    int DzubaRadialGrid(double r0=1.e-6, double rmax=250., double b=4.);
    int zeroNucleus();
    double enGuess(int Z, int n, int l, int tot_el, int num);





};

#endif
