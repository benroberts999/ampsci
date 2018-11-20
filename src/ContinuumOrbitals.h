#ifndef _CNTMORB_H
#define _CNTMORB_H
#include <string>
#include <vector>
#include <cmath>
#include "FPC_physicalConstants.h"
#include "ATI_atomInfo.h"
#include "ElectronOrbitals.h"
#include "ADAMS_solveLocalBS.h"
#include "ADAMS_solveLocalContinuum.h"


class ContinuumOrbitals{

  public:

    ContinuumOrbitals(ElectronOrbitals wf, int izion=1);
    //takes in grid, v from here

    int solveLocalContinuum(double ec, int min_l, int max_l);
    int solveLocalContinuum(double ec, int max_l);

    int solveZeffContinuum(double ec, double Zeff, int min_l, int max_l);

    void clear();

    std::vector< std::vector<double> > p;
    std::vector< std::vector<double> > q;
    std::vector<double> en;
    std::vector<int> kappa;

  private:

    std::vector<double> v; // v_nuc + v_dir

    int Z,Zion; //XXX Zion always 1 for now???

    double alpha;

    std::vector<double> r;
    std::vector<double> drdt;
    std::vector<double> dror;
    int NGPb;
    double h;


};

#endif
