#pragma once
#include "HF/YkTable.hpp" //for m_Yab
#include <string>
#include <vector>
class Wavefunction;
class DiracSpinor;
class Grid;

// Print-outs (for debugging) - work better without OMP
constexpr bool print_final_eps = false;
constexpr bool print_each_eps = false;

struct EpsIts {
  double eps;
  int its;
};

enum class HFMethod { HartreeFock, ApproxHF, Hartree, GreenPRM, TietzPRM };

//******************************************************************************
// For non-constant damping
// Slowly ramps the damping factor from a_beg to a_end over interval (beg, end)
static inline auto rampedDamp(double a_beg, double a_end, int beg, int end) {
  return [=](int i) {
    if (i >= end)
      return a_end;
    if (i <= beg)
      return a_beg;
    return (a_end * (i - beg) + a_beg * (end - i)) / (end - beg);
  };
}

//******************************************************************************

//! @brief
//! Solves relativistic Hartree-Fock equations

/*! @details
\par Construction
Requires Wavefunction object, and string core configuration (e.g., '[Xe],6s2')
\par Usage
Will solve HF equations for the core when constructed.
Note: it constructs the wf.core-orbitals vector... this is probably not the
right way to do it.

\par Definitions:

v^k_ab(r)   := Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
rho(r')     := fa(r')*fb(r') + ga(r')gb(r')
Lambda^k_ab := 3js((ja,jb,k),(-1/2,1/2,0))^2 * parity(la+lb+k)
vex[a]      := [v_ex*Fa](r) *(Fa/Fa^2) (approx exchange)
*/
class HartreeFock {

public:
  //! @brief HFMethod is enum class, eps_HF si convergence goal.
  //! @details h_d & g_t are parametric potential parameters (only used if
  //! method=Green/Teitz). Note: solves for HF core obitals (+populates them) on
  //! construct
  HartreeFock(HFMethod method, Wavefunction &wf, const std::string &in_core,
              double eps_HF = 0.0, double h_d = 0.0, double g_t = 0.0);

  //! @brief Solves HF for valence list; valence states must already be present
  //! in WaveFunction (bad, instead, give it a vector of DiracSpinors!)
  void solveValence();

  double calculateCoreEnergy() const;

public:
  //! @brief Calculates V_exch * Fa, for any orbital Fa (calculated Coulomb
  //! integral).
  //! @details Needs existing orbital Fa, and the core orbitals.
  //! k_cut is max multipolarity to sum over for exchange term [can limit to ~1
  //! (e.g.) for speed when high accuracy is not required]
  static DiracSpinor vex_psia_any(const DiracSpinor &Fa,
                                  const std::vector<DiracSpinor> &core,
                                  int k_cut = 99);

  bool verbose = true;

private:
  Wavefunction *const p_wf;
  const Grid *const p_rgrid;
  Coulomb::YkTable m_Yab;

  static constexpr int m_max_hf_its = 99;

  // Optionally force orthogonalisation. False by default.
  static constexpr bool m_explicitOrthog_cc = false;
  static constexpr bool m_explicitOrthog_cv = false;

  const double m_eps_HF;
  const bool m_excludeExchange;
  const HFMethod m_method;

public:
  bool excludeExchangeQ() const { return m_excludeExchange; }

private:
  // The "localised"/approximate HF potential:
  std::vector<std::vector<double>> appr_vex_core = {};
  std::vector<std::vector<double>> appr_vex_val = {};

private:
  void hf_core_approx(const double eps_target_HF);
  void starting_approx_core(const std::string &in_core, int log_converge = 3,
                            HFMethod method = HFMethod::GreenPRM,
                            double h_g = 0, double d_t = 0);

  DiracSpinor vex_psia(const DiracSpinor &Fa) const;
  void vex_psia(const DiracSpinor &Fa, DiracSpinor &vexFa) const;
  const std::vector<double> &get_vex(const DiracSpinor &psi) const;

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

public:
  HartreeFock &operator=(const HartreeFock &) = delete; // copy assignment
  HartreeFock(const HartreeFock &) = default;           // copy constructor
  ~HartreeFock() = default;

public:
  //! @brief Convers string (name) of method (HartreeFock, Hartree etc.) to enum
  static HFMethod parseMethod(const std::string &in_method);

  static std::vector<double>
  form_approx_vex_any(const DiracSpinor &Fa,
                      const std::vector<DiracSpinor> &core, int k_cut = 99);
};
