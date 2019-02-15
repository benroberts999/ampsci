#pragma once
#include "Grid.h"
#include <string>
#include <vector>

// enum class Operator { unity, gamma0, gamma5, dr, dr2 };
// enum class NucleusType { Fermi, spherical, zero };

struct Orbital {

  Orbital(int n, int k);

  std::vector<double> f;
  std::vector<double> g;
  double en;

  const int n;
  const int k;
  // const int state_index;

  double oc_frac;

  int its;
  double eps;
};
