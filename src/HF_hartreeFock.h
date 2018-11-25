#pragma once
#include <vector>
#include "ElectronOrbitals.h"

namespace HF{

  const int MAX_HART_ITS=100; //Max number of Hartree iterations
  const double default_eps=1.e-6;

  int hartreeFockCore(ElectronOrbitals &wf, double eps_HF);
  int hartreeFockValence(ElectronOrbitals &wf, int na, int ka, double eps_HF);
  int hartreeCore(ElectronOrbitals &wf, double eps_hartree=default_eps);

  int formNewVdir(const ElectronOrbitals &wf, std::vector<double> &vdir_new,
    bool core);

  int formVexCore(const ElectronOrbitals &wf,
    std::vector< std::vector<double> > &vex);
  int formVexA(const ElectronOrbitals &wf, int a, std::vector<double> &vex_a);
  int formLambdaABk(std::vector<double> &L_abk, int tja, int tjb, int la,
    int lb);

}
