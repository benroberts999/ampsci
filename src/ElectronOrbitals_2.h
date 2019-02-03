#pragma once
#include <string>
#include <vector>

enum class Operator { unity, gamma0, gamma5, dr, dr2 };
enum class NucleusType { Fermi, spherical, zero };

class Orbital {

  double *f;
  double *g; //???
  // also: max, pinf? (f, or f+g?)

  double en;

  const int n;
  const int k;
  const int state_index;

  double oc_frac;

  int its;
  double eps;
};

//******************************************************************************
class ElectronOrbitals {

private:
  std::vector<double> f_block;
  std::vector<double> g_block;

  Grid grid;

public:
  std::vector<Orbital> phi;
};
