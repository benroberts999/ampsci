#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>
class Wavefunction;
class Grid;

class ContinuumOrbitals {

public:
  ContinuumOrbitals(const Wavefunction &wf, int izion = 1);
  ContinuumOrbitals &operator=(const ContinuumOrbitals &) = delete;
  ContinuumOrbitals(const ContinuumOrbitals &) = default;
  ~ContinuumOrbitals() = default;
  // takes in grid, v from here

  int solveLocalContinuum(double ec, int min_l, int max_l);
  int solveLocalContinuum(double ec, int max_l);

  void clear();

  std::vector<DiracSpinor> orbitals = {};

private:
  const Grid *const p_rgrid;

  const int Z;
  const int Zion; // XXX Zion always 1 for now???
  const double alpha;

  std::vector<double> v = {}; // v_nuc + v_dir
};
