#pragma once
#include "Coulomb/YkTable.hpp"
#include "HF/Breit.hpp"
#include "Physics/Parametric_potentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/RadPot.hpp"
#include <memory>
#include <string>
#include <vector>
class Wavefunction;
class DiracSpinor;
class Grid;
namespace MBPT {
class CorrelationPotential;
}

//! Functions and classes for Hartree-Fock
namespace HF {

//==============================================================================
//! Small struct to store: {eps, its, symbol}. eps=convergence; its=iterations;
//! symbol=which state. May be sorted (by eps).
struct EpsIts {
  double eps{0.0};
  int its{0};
  std::string symbol{};
  friend bool operator<(const EpsIts &l, const EpsIts &r) {
    return l.eps < r.eps;
  }
};

//==============================================================================
//! @brief Methods available for self-consistant field model
/*! @details
 - HartreeFock: Self-consistent Hartree-Fock method
 - ApproxHF   : Approximate (localised) Hartree-Fock method
 - Hartree    : Core-Hartree method. No exchange, Vdir includes self-interaction
 - KohnSham   : Kohn-Sham (Density functional), includes Latter correction
 - Local      : Uses a local parameteric potential. [NOT self-consistant field]
 */
enum class Method { HartreeFock, ApproxHF, Hartree, KohnSham, Local };
//! @brief Convers string (name) of method (HartreeFock, Hartree etc.) to enum
Method parseMethod(const std::string &in_method);
std::string parseMethod(const Method &in_method);

//==============================================================================
//! Forms approx (localised) exchange potential, from scratch
//! @details Needs existing orbital Fa, and the core orbitals.
//! k_cut is max multipolarity to sum over for exchange term [can limit to ~1
//! (e.g.) for speed when high accuracy is not required]
std::vector<double> vex_approx(const DiracSpinor &Fa,
                               const std::vector<DiracSpinor> &core,
                               int k_cut = 99, double lambda_cut = 0.003);

//! @brief Calculates V_exch * Fa, for any orbital Fa (calculates Coulomb
//! integral from scratch).
//! @details  k_cut is max multipolarity to sum over for exchange term [can
//! limit to ~1 (e.g.) for speed when high accuracy is not required]
DiracSpinor vexFa(const DiracSpinor &Fa, const std::vector<DiracSpinor> &core,
                  int k_cut = 99);

//==============================================================================
//==============================================================================
//==============================================================================

//! Solves relativistic Hartree-Fock equations for core and valence. Optionally
//! includes Breit and QED effects. Can include Sigma (correlations) for valence
//! states. Class stores nuc. and direct potentials, a set of yk integrals, and
//! QED potential. Stores the core orbitals.
class HartreeFock {

private:
  std::shared_ptr<const Grid> m_rgrid;
  std::vector<DiracSpinor> m_core;
  std::vector<double> m_vnuc;
  std::optional<QED::RadPot> m_vrad;
  std::optional<HF::Breit> m_VBr;
  double m_alpha;
  Method m_method;
  double m_eps_HF;
  std::vector<double> m_vdir;
  Coulomb::YkTable m_Yab;
  int m_max_hf_its = 128;

public:
  //! @brief Method is enum class, eps_HF is convergence goal.
  /*! @details
    Required:
      - rgrid: Radial grid (shared pointer)
        - (This is required to allow no core orbitals)
        - Assumed to be same grid as for core orbitals
      - vnuc - nuclear potential
        - Assumed to be same length as radial grid
        - A copy is stored. May be updated
    Optional:
      - alpha (fine structure constant). default = true value
      - method default = HartreeFock
      - x_Breit - Breit scaling factor. 0=no Breit (default), 1=Breit. may set
        small number to check for non-linear contributions
      - eps_HF: convergence goal
      - potential: which parametric potential used for initial Potential
      - h and d (or g and t) are parameters for above (if left zero, default
        will be chosen)
      - Note:  Parametric::Type potential (and parameters H,d) are for
    initial approx. Usually doesn't matter at all, and defaults should be used.
    If using local potential [method=Local], these are the final parameters. If
    any are set to zero - will be looked up.
  */
  HartreeFock(std::shared_ptr<const Grid> rgrid, std::vector<double> vnuc,
              std::vector<DiracSpinor> core,
              std::optional<QED::RadPot> vrad = std::nullopt,
              double m_alpha = PhysConst::alpha,
              Method method = Method::HartreeFock, double x_Breit = 0.0,
              double eps_HF = 0.0,
              Parametric::Type potential = Parametric::Type::Green,
              double H_g = 0.0, double d_t = 0.0);

  //! Solves HF equations self-consitantly for core orbs. Returns epsilon.
  EpsIts solve_core(bool print = true);

  //! Solves HF for given valence list. They need not already be solutions.
  //! @details Note: If given energy is set to zero, states assumed to not be
  //! existing solutions; initial energy is guessed and solved from scratch. If
  //! initial energy is non-zero, that energy is used and states are assumed to
  //! already be (approximate) solutions.
  void
  solve_valence(std::vector<DiracSpinor> *valence, bool print = true,
                const MBPT::CorrelationPotential *const Sigma = nullptr) const;

  //! Solves HF equation (+ Sigma) for single valence state.
  EpsIts hf_valence(DiracSpinor &Fv,
                    const MBPT::CorrelationPotential *const Sigma = nullptr,
                    std::optional<double> eta = {},
                    std::optional<int> prev_its = {}) const;

  //! Calculates the HF core energy (not including Breit?)
  double calculateCoreEnergy() const;

  //! Calculates exchange term Vex*Fa
  DiracSpinor vexFa(const DiracSpinor &Fa) const {
    // calls static version with HF core
    return ::HF::vexFa(Fa, m_core, 99);
  }

  //! Breit interaction V_Br*Fa
  DiracSpinor VBr(const DiracSpinor &Fv) const;

  //---------------------------

  //! Resturns a const reference to the radial grid
  const Grid &grid() const { return *m_rgrid; };
  //! Resturns copy of shared_ptr to grid [shared resource] - used when we want
  //! to construct a new object that shares this grid
  std::shared_ptr<const Grid> grid_sptr() const { return m_rgrid; };

  //! Returns reference to Vdir (direct HF potential)
  const std::vector<double> &vdir() const { return m_vdir; }
  std::vector<double> &vdir() { return m_vdir; }

  //! Returns reference to Vnuc (nuclear potential)
  const std::vector<double> &vnuc() const { return m_vnuc; }
  std::vector<double> &vnuc() { return m_vnuc; }

  //! Electric part of radiative potential
  std::vector<double> Hrad_el(int l = 0) const;
  //! Magnetic (off-diagonal) part of radiative potential. Doesn't currently
  //! depend on l
  std::vector<double> Hmag(int l = 0) const;

  //! vlocal = vnuc + vrad_el + vdir
  std::vector<double> vlocal(int l = 0) const;

  //! Which method used to solve HF
  Method method() const { return m_method; }

  //! Effective charge at large Z : zion = Z - num_core_electrons
  double zion() const;

  //! Returns true if exchange not included
  bool excludeExchangeQ() const {
    return !(m_method == Method::HartreeFock || m_method == Method::ApproxHF);
  }

  //! vector of core orbitals
  const std::vector<DiracSpinor> &core() const { return m_core; }

  //! Value of fine-structure constant used
  double alpha() const { return m_alpha; }

  //! Update the Vrad used inside HF (only used if we want QED into valence but
  // not core, for testing)
  void set_Vrad(QED::RadPot in_vrad) { m_vrad = std::move(in_vrad); } // XXX
  //! Get (const) ptr to Vrad - may be null
  const QED::RadPot *Vrad() const { return m_vrad ? &*m_vrad : nullptr; }
  QED::RadPot *Vrad() { return m_vrad ? &*m_vrad : nullptr; }

  //! pointer to Breit - may be nullptr if no breit
  const HF::Breit *vBreit() const { return m_VBr ? &*m_VBr : nullptr; }
  //! Breit scale factor (usualy 0 or 1)
  double x_Breit() const { return m_VBr ? m_VBr->scale_factor() : 0.0; }

  //! Number of electrons in the core
  int num_core_electrons() const;

private:
  // Solve Dirac equation for core states (just once) using existing vdir
  // (usually, vdir set to parametric potential beforehand)
  EpsIts solve_initial_core(const double eps);
  // Solve equations self-consistantly for core, using local method (either
  // core-Hartree or Kohn-Sham)
  EpsIts selfcon_local_core(const double eps_target_HF);
  // Solve HF equations self-consistantly for core, using approximate HF method
  EpsIts hf_approx_core(const double eps_target_HF);
  // Solve HF equations self-consistantly for core, using Hartree-Fock method
  EpsIts hartree_fock_core();
  // Solves HF equation for valence state, assuming local potential (Hartree,
  // Local, Kohn-Sham or approxHF)
  EpsIts local_valence(DiracSpinor &Fa) const;

  /*
    // same as hf_valence, but uses Green method
    EpsIts hf_valence_Green(
        DiracSpinor &Fv,
        const MBPT::CorrelationPotential *const Sigma = nullptr) const;
  */

  // Solves HF equation for given state, using non-local Green's method for
  // inhomogeneous ODE (used for hartree_fock_core()).
  // Solve Dirac Equation (Eigenvalue):
  //  (H0 + Vl + Vx)Fa = 0
  //  (H0 + Vl)Fa = -VxFa
  // Vl is local (e.g., Vnuc + fVdir), Vx is non-local (e.g., (1-f)Vdir + Vex)
  // where v0 = (1-f)Vdir  [f=1 for valence states!, so v0 may be empty]
  // Vx also includes Breit, and Sigma
  // Small energy adjustmenets (and wfs), solve:
  // (Hl - e) dF = de * F -VxFa
  // e -> e+de, F->F+dF
  // Core is input so can call in a thread-safe way! (with a 'old_core' copy)
  // Only used in dE from dF
  void hf_orbital_green(
      DiracSpinor &Fa, double en, const std::vector<double> &vl,
      const std::vector<double> &H_mag, const DiracSpinor &VxF,
      const std::vector<DiracSpinor> &static_core,
      const std::vector<double> &dv0 = {}, const HF::Breit *const VBr = nullptr,
      const MBPT::CorrelationPotential *const Sigma = nullptr) const;

  // Calc's Vex*Fa, for Fa in the core. Fa must be in the core
  void vex_Fa_core(const DiracSpinor &Fa, DiracSpinor &vexFa) const;

  // Option to re-scale diract potential so that V(r)~-zion/r at large r
  enum class ReScale { yes = true, no = false };
  // Forms direct potential
  void update_vdir(ReScale re_scale = ReScale::no);
  // Adds the additional Kohn-Sham parts to Vdir
  void add_KohnSham_vdir_addition();

  // Sets Vdir to be parametric potential. By default, Greens potential
  void
  set_parametric_potential(bool print = true,
                           Parametric::Type potential = Parametric::Type::Green,
                           double H_g = 0.0, double d_t = 0.0);

  // Forms approximate vex for all core states
  void form_approx_vex_core(std::vector<std::vector<double>> &vex) const;
  std::vector<std::vector<double>> form_approx_vex_core() const;
  // Forms approximate vex for given core states
  void form_approx_vex_core_a(const DiracSpinor &Fa,
                              std::vector<double> &vex_a) const;
  std::vector<double> form_approx_vex_core_a(const DiracSpinor &Fa) const;

  // Energy guess for core states
  double enGuessCore(int n, int ka) const;
  // Energy guess for valence states
  double enGuessVal(int n, int ka) const;
};

} // namespace HF
