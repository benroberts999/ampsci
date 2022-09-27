#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <vector>
class Wavefunction;
class Grid;
namespace HF {
class HartreeFock;
}

//! Class stores set of continuum orbitals, and solves using Hartree-Fock method
class ContinuumOrbitals {

public:
  //! Takes existing Wavefunction object to construct. izion is effective charge at large r, izion = 1 usually (rarely matters).
  ContinuumOrbitals(const Wavefunction &wf, int izion = 1);
  //! Construct with Hartree-Fock pointer
  ContinuumOrbitals(const HF::HartreeFock *hf, int izion = 1);

  ContinuumOrbitals &operator=(const ContinuumOrbitals &) = default;
  ContinuumOrbitals(const ContinuumOrbitals &) = default;
  ~ContinuumOrbitals() = default;

  //! Solves continuum states with energy ec between min/max l
  int solveContinuumHF(double ec, int min_l, int max_l,
                       const DiracSpinor *psi = nullptr,
                       bool force_rescale = false, bool subtract_self = true,
                       bool force_orthog = true);

  //! Solves cntm states using simple H-like potential eith effective charge (Z_eff). Usually, Zeff = sqrt{2 * I_{njl}} * n.
  int solveContinuumZeff(double ec, int min_l, int max_l, double Z_eff,
                         const DiracSpinor *Fi, bool force_orthog);

  //! Checks orthogonality between cntm and core orbitals, returns worst eps
  double check_orthog(bool print = true) const;

  //! Resets (deletes) all orbitals
  void clear();

  std::vector<DiracSpinor> orbitals{};

private:
  void IncludeExchange(DiracSpinor &Fe, const DiracSpinor *psi,
                       bool force_orthog, const Grid &cgrid,
                       const std::vector<double> &vc, double r_asym);

  std::shared_ptr<const Grid> rgrid;
  const HF::HartreeFock *p_hf;
  int Zion;
  double alpha;
};
