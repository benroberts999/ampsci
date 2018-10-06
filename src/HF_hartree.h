#ifndef _HART_H
#define _HART_H
#include<cmath>
#include <string>
#include <vector>
#include "ElectronOrbitals.h"
#include "PRM_parametricPotentials.h"
#include "ATI_atomInfo.h"
#include "WIG_369j.h"
#include <fstream>

namespace HF{

  const int max_hartree=100; //Max number of Hartree iterations
  const double default_eps=1.e-6;

  int hartreeCore(ElectronOrbitals &wf, double eps_hartree=default_eps);
  
  int formNewVdir(ElectronOrbitals wf, std::vector<double> &vdir_new
  ,   bool core=true);

int formVexCore(ElectronOrbitals &wf, std::vector< std::vector<double> > &vex);
int formVexA(ElectronOrbitals &wf, int a, std::vector<double> &vex_a);
int formLambdaABk(std::vector<double> &L_abk, int tja, int tjb, int la, int lb);
double vexABr(ElectronOrbitals &wf, int a, int b, int ir,
  std::vector<double> &L_abk, int k_min);

}

#endif
