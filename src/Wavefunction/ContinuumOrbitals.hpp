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

  /*!
    @brief Solves continuum states with energy ec between min/max l
    @details
    If @p psi is given and @p subtract_self is true, the orbitals are solved
    in the V^{N-1} potential of the residual ion (hole in @p psi): the
    one-electron self-interaction V^psi_0 = D_psi - X_psi is removed from the
    Hartree-Fock potential. The direct part, y^0_{psi,psi}(r), is subtracted
    from the local potential (giving the Z_ion = 1 Coulomb tail); the exchange
    part (the one-electron self-exchange, @ref HF::vexFa_1el) is subtracted
    from the iterated non-local exchange.
  */
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
                       bool force_orthog, const std::vector<double> &vc,
                       bool subtract_self_exch = false);

  std::shared_ptr<const Grid> p_rgrid;
  const HF::HartreeFock *p_hf;
  double m_alpha;
};
