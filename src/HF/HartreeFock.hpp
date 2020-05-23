#pragma once
#include "Coulomb/YkTable.hpp" //for m_Yab
#include "Physics/PhysConst_constants.hpp"
#include <string>
#include <vector>
class Wavefunction;
class DiracSpinor;
class Grid;
namespace MBPT {
class CorrelationPotential;
}
namespace RadiativePotential {
class Vrad;
}

//! Functions and classes for Hartree-Fock
namespace HF {

// Print-outs (for debugging) - work better without OMP
constexpr bool print_final_eps = false;
constexpr bool print_each_eps = false;

//******************************************************************************
struct EpsIts {
  double eps;
  int its;
};

//******************************************************************************
enum class Method { HartreeFock, ApproxHF, Hartree };
//! @brief Convers string (name) of method (HartreeFock, Hartree etc.) to enum
Method parseMethod(const std::string &in_method);

//******************************************************************************
//! Forms approx (localised) exchange potential, from scratch
//! @details Needs existing orbital Fa, and the core orbitals.
//! k_cut is max multipolarity to sum over for exchange term [can limit to ~1
//! (e.g.) for speed when high accuracy is not required]
std::vector<double> vex_approx(const DiracSpinor &Fa,
                               const std::vector<DiracSpinor> &core,
                               int k_cut = 99);

//! @brief Calculates V_exch * Fa, for any orbital Fa (calculates Coulomb
//! integral from scratch).
//! @details Needs existing orbital Fa, and the core orbitals.
//! k_cut is max multipolarity to sum over for exchange term [can limit to ~1
//! (e.g.) for speed when high accuracy is not required]
DiracSpinor vexFa(const DiracSpinor &Fa, const std::vector<DiracSpinor> &core,
                  int k_cut = 99);

//******************************************************************************
//! @brief For non-constant damping. Slowly ramps the damping factor from a_beg
//! to a_end over interval (beg, end)
inline auto rampedDamp(double a_beg, double a_end, int beg, int end) {
  return [=](int i) {
    if (i >= end)
      return a_end;
    if (i <= beg)
      return a_beg;
    return (a_end * (i - beg) + a_beg * (end - i)) / (end - beg);
  };
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//! @brief Solves relativistic Hartree-Fock equations
/*! @details
\par Construction
Requires set of core orbitals + energies that are already solved! (it's the
energies that are most important, the energies must already be good guesses).
WaveFunction gives routine to solve initial core approx.

\par Usage
Solves HF equations for the core when solveCore() is called.
Core orbitals must already exist. Note: stores a pointer to external core
orbitals. These must not be extended/deleted while HF object exists.

\par Definitions:
v^k_ab(r)   := Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
rho(r')     := fa(r')*fb(r') + ga(r')gb(r')
Lambda^k_ab := 3js((ja,jb,k),(-1/2,1/2,0))^2 * parity(la+lb+k)
vex[a]      := [v_ex*Fa](r) *(Fa/Fa^2) (approx exchange)
*/
class HartreeFock {

public:
  //! @brief Method is enum class, eps_HF is convergence goal.
  HartreeFock(const Grid &in_grid, const std::vector<double> &in_vnuc,
              std::vector<DiracSpinor> *in_core,
              const RadiativePotential::Vrad *const in_vrad = nullptr,
              double m_alpha = PhysConst::alpha,
              Method method = Method::HartreeFock, double eps_HF = 0.0);
  HartreeFock(Wavefunction *wf, Method method = Method::HartreeFock,
              double eps_HF = 0.0);

  //! Solves HF equations self-consitantly for core orbs. Produces Vdir
  const std::vector<double> &solveCore();

  //! @brief Solves HF for valence list; valence states must already be present
  //! in WaveFunction (bad, instead, give it a vector of DiracSpinors!)
  void solveValence(std::vector<DiracSpinor> *valence, const bool print = true);

  //! Solves HF+Sigma equation: valence Brueckner orbitals. Writes to valence
  void solveBrueckner(std::vector<DiracSpinor> *valence,
                      const MBPT::CorrelationPotential &Sigma2,
                      const bool print = true);

  //! Calculates the HF core energy
  double calculateCoreEnergy() const;

  //! Returns const ref to V_dir (Direct HF potential)
  const std::vector<double> &get_vdir() const { return m_vdir; }

  //! Calculates exchange term Vex*Fa (calls static version with HF core)
  DiracSpinor calc_vexFa(const DiracSpinor &Fa) const {
    return vexFa(Fa, *p_core, 99);
  }
  //! Same as 'calc_vexFa'; Fa must be in core [y_ab exist already, so faster]
  DiracSpinor calc_vexFa_core(const DiracSpinor &Fa) const;

  //! Returns true of exchange not included
  bool excludeExchangeQ() const { return m_excludeExchange; }

  std::vector<double> get_vlocal(int l) const;

  int num_core_electrons() const;
  double get_alpha() const { return m_alpha; }
  const std::vector<DiracSpinor> &get_core() const { return *p_core; }

public:
  bool verbose = true; // update to input??
  const Grid *const p_rgrid;

private:
  const std::vector<double> *const p_vnuc;
  const RadiativePotential::Vrad *const p_vrad;

public:
  const double m_alpha;
  const Method m_method;

private:
  const double m_eps_HF;
  // pointer to core orbs that exist outside.
  // Note: core.size() must not change; core elements must remain in orig. order
  std::vector<DiracSpinor> *const p_core;
  std::vector<double> m_vdir;
  Coulomb::YkTable m_Yab;
  const bool m_excludeExchange; // XXX Kill this. Only HF,H,aHF

  static constexpr int m_max_hf_its = 99;
  // Optionally force orthogonalisation. False by default.
  static constexpr bool m_explicitOrthog_cc = false;
  static constexpr bool m_explicitOrthog_cv = false;

  const std::vector<double> &get_Hrad_el(int l) const;

public:
  const std::vector<double> &get_Hrad_mag(int l) const;

private:
  void hf_core_approx(const double eps_target_HF);
  void hf_core_refine();

  EpsIts hf_valence_approx(DiracSpinor &phi, double eps_target);
  EpsIts hf_valence_refine(DiracSpinor &phi);
  EpsIts hf_Brueckner(DiracSpinor &Fa, const MBPT::CorrelationPotential &Sigma);

  void hf_orbital(DiracSpinor &phi, double en, const std::vector<double> &vl,
                  const std::vector<double> &H_mag, const DiracSpinor &vx_phi,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<double> &v0 = {}) const;

  void brueckner_orbital(DiracSpinor &Fa, double en,
                         const std::vector<double> &vl,
                         const std::vector<double> &H_mag,
                         const DiracSpinor &VxF,
                         const MBPT::CorrelationPotential &Sigma,
                         const std::vector<DiracSpinor> &static_core) const;

  // Calc's Vex*Fa, for Fa in the core
  void vex_psia_core(const DiracSpinor &Fa, DiracSpinor &vexFa) const;
  // Forms direct potential
  void form_vdir(std::vector<double> &vdir, bool re_scale = false) const;
  // Forms approximate vex for all core states
  void form_approx_vex_core(std::vector<std::vector<double>> &vex) const;
  // Forms approximate vex for given core states
  void form_approx_vex_core_a(const DiracSpinor &Fa,
                              std::vector<double> &vex_a) const;
  std::vector<double> form_approx_vex_core_a(const DiracSpinor &Fa) const;

public:
  HartreeFock &operator=(const HartreeFock &) = delete; // copy assignment
  HartreeFock(const HartreeFock &) = default;           // copy constructor
  ~HartreeFock() = default;
};

} // namespace HF
