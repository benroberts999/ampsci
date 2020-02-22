#pragma once
#include "HF/CoulombIntegrals.hpp" //for m_cint
#include <string>
#include <vector>
class Wavefunction;
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
  vex[a]      := [v_ex*Fa](r) *(Fa/Fa^2) (approx exchange)
*/

struct EpsIts {
  double eps;
  int its;
};

enum class HFMethod { HartreeFock, ApproxHF, Hartree, GreenPRM, TietzPRM };

//******************************************************************************
class HartreeFock {
  friend class Coulomb;

public:
  static DiracSpinor
  solveMixedState(const DiracSpinor &Fa, const int k, const double omega,
                  const std::vector<double> &vl, const double alpha,
                  const std::vector<DiracSpinor> &core, const DiracSpinor &hFa);

public:
  static HFMethod parseMethod(const std::string &in_method);

  HartreeFock(HFMethod method, Wavefunction &wf, const std::string &in_core,
              double eps_HF = 0, double h_d = 0, double g_t = 0);

  // for HF basis:
  HartreeFock(Wavefunction &wf, const std::vector<DiracSpinor> &val_orbitals,
              double eps_HF = 0.0, bool in_ExcludeExchange = false);

  void solveValence();

  double calculateCoreEnergy() const;

  const std::vector<double> &get_vex(const DiracSpinor &psi) const;

  DiracSpinor vex_psia(const DiracSpinor &Fa) const;
  void vex_psia(const DiracSpinor &Fa, DiracSpinor &vexPsi) const;

  static DiracSpinor vex_psia_any(const DiracSpinor &Fa,
                                  const std::vector<DiracSpinor> &core,
                                  int k_cut = 99);
  static std::vector<double>
  form_approx_vex_any(const DiracSpinor &Fa,
                      const std::vector<DiracSpinor> &core, int k_cut = 99);

  bool verbose = true;

private:
  Wavefunction *const p_wf;
  const Grid *const p_rgrid;
  Coulomb m_cint;

  static constexpr int m_max_hf_its = 99;

  // Optionally force orthogonalisation. False by default.
  static constexpr bool m_explicitOrthog_cc = false;
  static constexpr bool m_explicitOrthog_cv = false;

  const double m_eps_HF;

public: // blahhhhh
  const bool m_excludeExchange;
  const HFMethod m_method;

private:
  // The "localised"/approximate HF potential:
  std::vector<std::vector<double>> appr_vex_core;
  std::vector<std::vector<double>> appr_vex_val;

private:
  void hf_core_approx(const double eps_target_HF);
  void starting_approx_core(const std::string &in_core, int log_converge = 3,
                            HFMethod method = HFMethod::GreenPRM,
                            double h_g = 0, double d_t = 0);

  void form_vdir(std::vector<double> &vdir, bool re_scale = false) const;
  void form_approx_vex_core(std::vector<std::vector<double>> &vex) const;
  void form_approx_vex_a(const DiracSpinor &Fa,
                         std::vector<double> &vex_a) const;

  void hf_core_refine();

  EpsIts hf_valence(DiracSpinor &phi, std::vector<double> &vexa);
  EpsIts hf_valence_approx(DiracSpinor &phi, std::vector<double> &vexa,
                           double eps_target_HF);
  EpsIts hf_valence_refine(DiracSpinor &phi);

  void hf_orbital(DiracSpinor &phi, double en, const std::vector<double> &vl,
                  const std::vector<double> &H_mag, const DiracSpinor &vx_phi,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<double> &v0 = {}) const;
};
