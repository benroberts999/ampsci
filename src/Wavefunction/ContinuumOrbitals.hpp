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
  //! Takes existing Wavefunction object to construct.
  ContinuumOrbitals(const Wavefunction &wf);
  //! Construct with Hartree-Fock pointer
  ContinuumOrbitals(const HF::HartreeFock *hf);

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
  void IncludeExchange(DiracSpinor &F_cntm, const DiracSpinor *F_i,
                       bool force_orthog, const std::vector<double> &vc);

  std::shared_ptr<const Grid> p_rgrid;
  const HF::HartreeFock *p_hf;
  double m_alpha;
};
