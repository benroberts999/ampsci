#ifndef _HART_H
#define _HART_H
#include<cmath>
#include <string>
#include <vector>
#include "ElectronOrbitals.h"
#include "PRM_parametricPotentials.h"
#include <fstream>

const int max_hartree=100; //Max number of Hartree iterations
const double default_eps=1.e-6;

int HF_hartreeCore(ElectronOrbitals &wf, double eps_hartree=default_eps);
int formNewVdir(ElectronOrbitals wf, std::vector<double> &vdir_new, bool core=true);


#endif
