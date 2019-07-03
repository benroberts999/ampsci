#pragma once
#include "CoulombIntegrals.hpp"
#include <string>
#include <vector>
class ElectronOrbitals;
class DiracSpinor;
class Grid;

/*
Calculates self-consistent Hartree-Fock potential, including exchange.
Solves all core and valence states.

Note: can run without exchange (Hartree) - for tests only.
(Still, non-local potential, 'vex' different for each core state)

XXX Have option to give a list of valence states!
Can solve them to some degree in parallel
Requires re-writing the valence part (a little)

  Definitions:
  v^k_ab(r)   := Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
  rho(r')     := fa(r')*fb(r') + ga(r')gb(r')
  Lambda^k_ab := 3js((ja,jb,k),(-1/2,1/2,0))^2 * parity(la+lb+k)
  vex[a]      := [v_ex*psi_a](r) *(psi_a/psi_a^2) (approx exchange)
*/

class HartreeFock {
  friend class Coulomb;

public:
  HartreeFock(ElectronOrbitals &wf, const std::string &in_core,
              double eps_HF = 1.e-8, bool in_ExcludeExchange = false);

  void solveNewValence(int n, int kappa);
  void solveValence(DiracSpinor &phi, std::vector<double> &vexa);

  double calculateCoreEnergy() const;

  const std::vector<double> &get_vex(const DiracSpinor &psi) const;

  bool verbose = true;

private:
  ElectronOrbitals *const p_wf;
  const Grid *const p_rgrid;

  Coulomb m_cint;

  double m_eps_HF = 1.e-8;

  static const int MAX_HART_ITS = 99;
  const bool m_excludeExchange; // for testing

  // The "localised"/approximate HF potential:
  std::vector<std::vector<double>> vex_core;
  std::vector<std::vector<double>> vex_val;

  std::vector<std::vector<std::vector<double>>> m_arr_Lambda_nmk;

private:
  void hartree_fock_core();
  void starting_approx_core(const std::string &in_core);

  void form_vdir(std::vector<double> &vdir, bool re_scale = false) const;
  void form_approx_vex_core(std::vector<std::vector<double>> &vex) const;
  void form_approx_vex_a(const DiracSpinor &phi_a,
                         std::vector<double> &vex_a) const;
};
