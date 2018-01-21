#ifndef _ORBITALS_H
#define _ORBITALS_H
#include <string>
#include <vector>
#include <cmath>
#include "ElectronOrbitals.h"
#include "physicalConstants.h"
#include "atomInfo.h"
#include "adamsSolveLocalBS.h"

const int NGP_DEFAULT=1000; //???

class ElectronOrbitals{

  public:

    ElectronOrbitals(int in_z, int in_a=0, int in_ngp=NGP_DEFAULT);
    ElectronOrbitals(std::string s_in_z, int in_a=0, int in_ngp=NGP_DEFAULT);
    // z
    // ion=0, a, ngp, num_pot_type

    std::vector< std::vector<double> > f;
    std::vector< std::vector<double> > g;
    std::vector<double> en;

    std::vector<double> r;
    std::vector<double> drdt;
    std::vector<double> dror;

    std::vector<double> vnuc;

    int Z,A;
    std::string atom;

    double var_alpha; // like this?

    int max_n;
    int max_l;
    int num_states; //?

    std::vector<unsigned> nlist;
    std::vector<int> klist;
    //std::vector<unsigned int> llist;
    //std::vector<unsigned int> j2list;
    //std::vector<std::string> hrlist;

    localBoundState();


  private:

    //Number of grid points:
    int ngp;

    int formRadialGrid();
    int zeroNucleus();
    int sphericalNucleus(double rnuc=0);
    int fermiNucleus(double t=0, double c=0);


};

#endif
