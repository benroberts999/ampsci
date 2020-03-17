#pragma once
#include "HF/HartreeFockClass.hpp" // forward decl..
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp" // NonRelSEConfig
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp" //PhysConst::alpha
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static bool dummy_bool{};
//******************************************************************************
/*!
@brief Stores Wavefunction (set of core+valence orbitals, grid etc.)
@details
\par Construction:
  - Needs Z (atom); either as int or string
  - Set of GridParameters [see Maths/Grid]
  - Set of Nuclear::Parameters [see Physics/NuclearPotentials]
  - var_alpha = \f$\lambda\f$, \f$\alpha = \lambda\alpha_0\f$

Note: cannot be copied. Hope to remedy this soon; would be nice.
*/
class Wavefunction {

public:
  template <typename T>
  Wavefunction(T in_z, const GridParameters &gridparams,
               const Nuclear::Parameters &nuc_params, double var_alpha = 1.0)
      : rgrid({gridparams}),                                        //
        m_alpha(PhysConst::alpha * var_alpha),                      //
        m_Z(AtomData::get_z(in_z)),                                 //
        m_A(nuc_params.a),                                          //
        m_nuc_params(nuc_params),                                   //
        vnuc(Nuclear::formPotential(nuc_params, m_Z, m_A, rgrid.r)) //
  {
    if (m_alpha * m_Z > 1.0) {
      std::cerr << "Alpha too large: Z*alpha=" << m_Z * m_alpha << "\n";
      std::abort();
    }
  }

public:
  std::vector<DiracSpinor> core_orbitals = {};
  std::vector<DiracSpinor> valence_orbitals = {};
  std::vector<DiracSpinor> basis = {};
  const Grid rgrid;

private:
  const double m_alpha; // store internal value for alpha (allows variation)
  const int m_Z, m_A;
  Nuclear::Parameters m_nuc_params;
  std::unique_ptr<HF::HartreeFock> m_pHF = nullptr;

public:
  //! Nuclear potential
  std::vector<double> vnuc = {};
  //! Direct/local part of the electron potential
  std::vector<double> vdir = {}; //
  //! QED magnetic form-factor
  std::vector<double> Hse_mag = {}; // magnetic form-factor

private:
  // Core configuration (non-rel terms)
  std::vector<AtomData::NonRelSEConfig> m_core_configs = {};
  int num_core_electrons = 0; // Nc = N - M
  std::string m_core_string = "";

public: // const methods: "views" into WF object
  // Rule is: if function is single-line, define here. Else, in .cpp
  double get_alpha() const { return m_alpha; }
  int Znuc() const { return m_Z; }
  int Anuc() const { return m_A; }
  //! Number of neutrons, A-Z
  int Nnuc() const { return (m_A > m_Z) ? (m_A - m_Z) : 0; }
  int Ncore() const { return num_core_electrons; }
  const Nuclear::Parameters &get_nuclearParameters() const {
    return m_nuc_params;
  }
  bool exclude_exchangeQ() const {
    if (m_pHF == nullptr)
      return true;
    return m_pHF->excludeExchangeQ();
  }

  std::size_t getStateIndex(int n, int k, bool &is_valence = dummy_bool) const;

  //! Finds requested state; returns nullptr if not found
  //! @details is_valence is optional out-parameter; tells you where orb was
  //! found
  const DiracSpinor *getState(int n, int k,
                              bool &is_valence = dummy_bool) const;

  //! Returns full core configuration
  std::string coreConfiguration() const { return m_core_string; }
  //! Returns core configuration, in nice output notation
  std::string coreConfiguration_nice() const {
    return AtomData::niceCoreOutput(m_core_string);
  }
  //! Outputs screen-friendly nuclear parameters
  std::string nuclearParams() const;

  //! Returns string of atom info (Z, A)
  std::string atom() const {
    return AtomData::atomicSymbol(m_Z) + ", Z=" + std::to_string(m_Z) +
           " A=" + std::to_string(m_A);
  }
  std::string atomicSymbol() const { return AtomData::atomicSymbol(m_Z); }

  //! Prints table of core orbitals + energies etc. Optionally sorted by energy
  void printCore(bool sorted = true) const;

  //! @breif Prints table of valence orbitals + energies etc. Optionally sorted
  //! by energy
  //! @details Can optionally give it any list of orbitals to print
  void printValence(bool sorted = true,
                    const std::vector<DiracSpinor> &tmp_orbitals = {}) const;
  //! Prints table of Basis orbitals, compares to HF orbitals
  void printBasis(bool sorted = false) const;
  bool isInCore(int n, int k) const;
  bool isInValence(int n, int k) const;
  bool isInCore(const DiracSpinor &phi) const;
  int maxCore_n(int ka_in = 0) const;
  int maxCore_l() const;

  //! Calculated rho(r) = sum_c psi^2(r) for core states
  std::vector<double> coreDensity() const;

  //! Performs hartree-Fock procedure for core: note: poplulates core
  void hartreeFockCore(HF::Method method, const std::string &in_core,
                       double eps_HF = 0, double h_d = 0, double g_t = 0);

  auto coreEnergyHF() const;

  //! Performs hartree-Fock procedure for valence: note: poplulates valnece
  void hartreeFockValence(const std::string &in_valence_str);

  //! Calculates radiative potential. Stores in vnuc, and Hmag
  void radiativePotential(double x_Ueh, double x_SEe_h, double x_SEe_l,
                          double x_SEm, double rcut, double scale_rN);

  //! Calculates + populates basis [see BSplineBasis]
  void formBasis(const std::string &states_str, const std::size_t n_spl,
                 const std::size_t k_spl, const double r0_spl,
                 const double rmax_spl, const bool positronQ = false);

  //! @brief Solves Dirac bound state problem, with optional 'extra' potential
  //! log_eps is log_10(convergence_target).
  void solveDirac(DiracSpinor &psi, double e_a, const std::vector<double> &vex,
                  int log_eps = 0) const;
  void solveDirac(DiracSpinor &psi, double e_a = 0, int log_eps = 0) const;

  //! Populates core orbitals accorind to given core string (+solves)
  void solveInitialCore(const std::string &str_core_in, int log_dele_or = 0);
  //! Adds new valence orbtial (+solves using vdir)
  void solveNewValence(int n, int k, double en_a = 0, int log_dele_or = 0);

  //! energy guess for a core state with n,l; quite rough
  double enGuessCore(int n, int l) const;
  //! energy guess for a valence state with n,l
  double enGuessVal(int n, int ka) const;

  //! (approximately) OrthoNormalises a set of any orbitals.
  static void orthonormaliseOrbitals(std::vector<DiracSpinor> &in_orbs,
                                     int num_its = 1);
  //! (exactly) OrthoNormalises psi_v against of any orbitals.
  static void orthonormaliseWrt(DiracSpinor &psi_v,
                                const std::vector<DiracSpinor> &in_orbs);
  //! (exactly) OrthoGonalises psi_v against of any orbitals (no Norm).
  static void orthogonaliseWrt(DiracSpinor &psi_v,
                               const std::vector<DiracSpinor> &in_orbs);

private:
  void determineCore(std::string str_core_in);
  // not static, since skip_core. Make static version??
  std::vector<AtomData::DiracSEnken>
  listOfStates_nk(int num_val, int la, int lb = 0, bool skip_core = true) const;
  static std::vector<std::size_t>
  sortedEnergyList(const std::vector<DiracSpinor> &tmp_orbs,
                   bool do_sort = false);
};
