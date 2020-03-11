#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>
class Wavefunction;
class Grid;

//! Class similar to Wavefunction. Stores set of continuum orbitals
class ContinuumOrbitals {

public:
  //! Takes existing Wavefunction object to construct. izion = 1 usually
  ContinuumOrbitals(const Wavefunction &wf, int izion = 1);
  ContinuumOrbitals &operator=(const ContinuumOrbitals &) = delete;
  ContinuumOrbitals(const ContinuumOrbitals &) = default;
  ~ContinuumOrbitals() = default;
  // takes in grid, v from here

  //! Solves continuum states with energy ec between min/max l
  int solveLocalContinuum(double ec, int min_l, int max_l);
  //! Solves continuum states with energy ec between l=0 and l=max_l
  int solveLocalContinuum(double ec, int max_l);

  //! Resets (deletes) all orbitals
  void clear();

  std::vector<DiracSpinor> orbitals = {};

private:
  const Grid *const p_rgrid;

  const int Z;
  const int Zion; // XXX Zion always 1 for now???
  const double alpha;

  std::vector<double> v = {}; // v_nuc + v_dir
};
