#pragma once
#include "DiracSpinor.hpp"
#include <vector>
class ElectronOrbitals;
class Grid;

class ContinuumOrbitals {

public:
  ContinuumOrbitals(const ElectronOrbitals &wf, int izion = 1);
  // takes in grid, v from here

  int solveLocalContinuum(double ec, int min_l, int max_l);
  int solveLocalContinuum(double ec, int max_l);

  // int solveZeffContinuum(double ec, double Zeff, int min_l, int max_l);

  void clear();

  // XXX Swap to Dirac orbital!

  std::vector<DiracSpinor> orbitals;

  std::vector<std::vector<double>> f;
  std::vector<std::vector<double>> g;
  std::vector<double> en;
  std::vector<int> kappa;

private:
  const Grid *const p_rgrid;

  const int Z;
  const int Zion; // XXX Zion always 1 for now???
  const double alpha;

  std::vector<double> v; // v_nuc + v_dir

  // Do I need to store these? Grid? Or these elongated ones?
  // std::vector<double> r;
  // std::vector<double> drdt;
  // std::vector<double> dror;
  // std::size_t NGPb;
  // double h;
};
