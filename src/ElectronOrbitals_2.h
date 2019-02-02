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
struct Grid {

  Grid();

  std::vector<double> r;
  std::vector<double> drdu;
  std::vector<double> drduor;
  const double r0;
  const double rmax;
  const int ngp;
  const double du;

private:
  void exponentialRadialGrid(int ngp_in, double r0, double rmax);
  void logLinearRadialGrid(int ngp_in, double r0, double rmax, double b = 4.);
  void logLinearRadialGrid(double in_h, double r0, double rmax, double b = 4.);
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
