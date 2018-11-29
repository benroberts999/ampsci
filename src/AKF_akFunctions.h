#pragma once
#include "ElectronOrbitals.h"
#include "ExponentialGrid.h"
#include <vector>
#include <string>

namespace AKF{

  double CLkk(int L, int ka, int kb);

  void writeToTextFile(
    std::string fname,
    const std::vector< std::vector< std::vector<float> > > &AK,
    const std::vector<std::string> &nklst,
    double qmin, double qmax,
    double demin, double demax);

  int akReadWrite(std::string fname, bool write,
    std::vector< std::vector< std::vector<float> > > &AK,
    std::vector<std::string> &nklst,
    double &qmin, double &qmax,
    double &dEmin, double &dEmax);

  int calculateK_nk(const ElectronOrbitals &wf, int nk, int max_L, double dE,
    std::vector< std::vector<std::vector<float> > > &jLqr_f,
    std::vector<float> &K_nk, double Zeff=-1);

  int calculateKpw_nk(const ElectronOrbitals &wf, int nk, double dE,
    std::vector< std::vector<float> > &jl_qr,
    std::vector<float> &K_nk);

  void sphericalBesselTable(
    std::vector< std::vector< std::vector<float> > > &jLqr_f,
    int max_L,
    //double qmin, double qmax, int qsteps,
    const ExpGrid &qgrid,
    const std::vector<double> &r);

}
