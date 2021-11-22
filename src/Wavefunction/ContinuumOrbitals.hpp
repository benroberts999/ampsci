#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <vector>
class Wavefunction;
class Grid;
namespace HF {
class HartreeFock;
}

//! Class similar to Wavefunction. Stores set of continuum orbitals
class ContinuumOrbitals {

public:
  //! Takes existing Wavefunction object to construct. izion = 1 usually
  ContinuumOrbitals(const Wavefunction &wf, int izion = 1);
  ContinuumOrbitals &operator=(const ContinuumOrbitals &) = delete;
  ContinuumOrbitals(const ContinuumOrbitals &) = default;
  ~ContinuumOrbitals() = default;

  //! Solves continuum states with energy ec between min/max l
  int solveContinuumHF(double ec, int min_l, int max_l,
                       const DiracSpinor *Fi = nullptr);
  //! Solves continuum states with energy ec between l=0 and l=max_l
  int solveContinuumHF(double ec, int max_l, const DiracSpinor *Fi = nullptr);

  double check_orthog(bool print = true) const;

  //! Resets (deletes) all orbitals
  void clear();

  std::vector<DiracSpinor> orbitals{};

private:
  std::shared_ptr<const Grid> rgrid;
  const HF::HartreeFock *const p_hf;
  const int Z;
  const int Zion;
  const double alpha;
  std::vector<double> v_local; // v_nuc + v_dir
};
