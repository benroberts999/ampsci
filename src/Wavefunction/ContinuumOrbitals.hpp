#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <vector>
class Wavefunction;
class Grid;
namespace HF {
class HartreeFock;
}

//! Stores a set of continuum orbitals and solves them using the Hartree-Fock method
class ContinuumOrbitals {

public:
  //! Construct from an existing Wavefunction
  ContinuumOrbitals(const Wavefunction &wf);
  //! Construct from a HartreeFock pointer
  ContinuumOrbitals(const HF::HartreeFock *hf);

  ContinuumOrbitals &operator=(const ContinuumOrbitals &) = default;
  ContinuumOrbitals(const ContinuumOrbitals &) = default;
  ~ContinuumOrbitals() = default;

  //! Solves continuum states with energy ec for all kappa with l in [min_l, max_l].
  //! Fi: ionised orbital (used for self-interaction subtraction and projection).
  //! subtract_self: subtracts the hole-particle direct+exchange potential V_hp.
  //! force_orthog: enforces orthogonality to core after solving.
  int solveContinuumHF(double ec, int min_l, int max_l,
                       const DiracSpinor *Fi = nullptr,
                       bool force_rescale = false, bool subtract_self = true,
                       bool force_orthog = true);

  //! Solves continuum states using a simple H-like potential with effective
  //! charge Z_eff. Typically Z_eff = sqrt(2*I_{njl}) * n (DarkARC convention).
  int solveContinuumZeff(double ec, int min_l, int max_l, double Z_eff,
                         const DiracSpinor *Fi, bool force_orthog);

  //! Returns the worst |<continuum|core>| overlap; optionally prints all.
  double check_orthog(bool print = true) const;

  //! Clears all stored orbitals
  void clear();

  std::vector<DiracSpinor> orbitals{};

  //! V_hp|Fa> = (1/(2j_i+1)) vex(Fa,Fi) + y^0_{ii}·Fa
  static DiracSpinor hp_apply_V0(const DiracSpinor &Fa, const DiracSpinor *Fi,
                                 const std::vector<double> &vdir_0);

  //! Precomputes V_hp_core[a] = V_hp|a> and Vhp_ba[b,a] = <b|V_hp|a> for same-kappa core pairs.
  //! Results written via output pointers V_hp_core and Vhp_ba.
  static void hp_precompute(const std::vector<DiracSpinor> &core, int kappa,
                            const DiracSpinor *Fi,
                            const std::vector<double> &vdir_0,
                            std::vector<DiracSpinor> *V_hp_core,
                            std::vector<double> *Vhp_ba);

  //! Returns the Hermitian projection correction:
  //!   Pc V_hp|Fc> + V_hp Pc|Fc> - Pc V_hp Pc|Fc>
  //! where Pc projects onto same-kappa core states.
  static DiracSpinor
  add_hp_projection(const DiracSpinor &Fc, const DiracSpinor *Fi,
                    const std::vector<double> &vdir_0,
                    const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &V_hp_core,
                    const std::vector<double> &Vhp_ba);

private:
  void IncludeExchange(DiracSpinor *F_cntm, const DiracSpinor *F_i,
                       const std::vector<double> &vc,
                       const std::vector<double> &vdir_0);

  std::shared_ptr<const Grid> p_rgrid;
  const HF::HartreeFock *p_hf;
  double m_alpha;
};
