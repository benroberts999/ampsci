#ifndef _ORBITALS_H
#define _ORBITALS_H
#include <string> //???
#include <vector>
#include <cmath>
#include "physicalConstants.h"
#include "atomInfo.h"
#include "adamsSolveLocalBS.h"
#include <gsl/gsl_sf_fermi_dirac.h>

const int NGP_DEFAULT=1000; //???

class ElectronOrbitals{

  public:

    ElectronOrbitals(int in_z, int in_a, int in_ngp=NGP_DEFAULT,
        double var_alpha=1);

    std::vector< std::vector<double> > p;
    std::vector< std::vector<double> > q;
    std::vector<double> en;

    std::vector<double> r;
    std::vector<double> drdt;
    std::vector<double> dror;
    int ngp;
    double h;

    std::vector<double> vnuc;

    int Z,A;
    std::string atom;

    //double var_alpha; // like this?
    double alpha;

    int max_n;
    int max_l;
    int num_states; //?

    std::vector<unsigned> nlist;
    std::vector<int> klist;

    std::vector<int> pinflist;
    std::vector<int> itslist;
    std::vector<double> epslist;

    int localBoundState(int in_max_n, int in_max_l=100);

    double diracen(int z, int n, int k);
    int sphericalNucleus(double rnuc=0);
    int fermiNucleus(double t=0, double c=0);


  private:

    //Number of grid points:


    int formRadialGrid();
    int zeroNucleus();





};

#endif
