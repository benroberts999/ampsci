#pragma once
// #include "Grid.h"
#include <string>
#include <vector>

// enum class Operator { unity, gamma0, gamma5, dr, dr2 };
// enum class NucleusType { Fermi, spherical, zero };

struct DiracSpinor {

  DiracSpinor(int in_n, int in_k, int ngp)
      : en(0), n(in_n), k(in_k), pinf(ngp), its(-1), eps(-1.), occ_frac(-1.)
  //
  {
    f.resize(ngp, 0);
    g.resize(ngp, 0);
  }

  std::vector<double> f;
  std::vector<double> g;
  double en;

  const int n;
  const int k;

  int pinf;
  int its;
  double eps;
  double occ_frac;
};
