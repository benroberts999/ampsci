#ifndef _AKFUNS_H
#define _AKFUNS_H
#include "ElectronOrbitals.h"
#include "ContinuumOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include "HF_hartree.h"
#include "WIG_369j.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <gsl/gsl_sf_bessel.h>


//Declare angular coeficient. [See Phys.Rev.D 93, 115037 (2016).]
double CLkk_OLD(int L, int ka, int kb);
double CLkk(int L, int ka, int kb);

int akReadWrite(std::string fname, bool write,
  std::vector< std::vector< std::vector<float> > > &AK,
  std::vector<float> &dElst,
  std::vector<std::string> &nklst,
  std::vector<float> &qlst);

int calculateK_nk(ElectronOrbitals &wf, int nk, int max_L, double dE,
  std::vector< std::vector<std::vector<float> > > &jLqr_f,
  std::vector< std::vector<float> > &K_nk);

int calculateKpw_nk(ElectronOrbitals &wf, int nk, double dE,
  std::vector< std::vector<float> > &jl_qr,
  std::vector< std::vector<float> > &K_nk
);

#endif
