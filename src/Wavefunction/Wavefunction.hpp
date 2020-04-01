#pragma once
#include "HF/HartreeFockClass.hpp" // forward decl..
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp" // NonRelSEConfig
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp" //PhysConst::alpha
#include "Physics/RadiativePotential.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <utility>
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
  std::vector<DiracSpinor> core_orbitals{};
  std::vector<DiracSpinor> valence_orbitals{};
  //! Basis, daigonalised over HF core. Used for MBPT
  std::vector<DiracSpinor> basis{};
  //! Sprectrum: like basis, but includes Sigma.
  std::vector<DiracSpinor> spectrum{};
  const Grid rgrid;

private:
  const double m_alpha; // store internal value for alpha (allows variation)
  const int m_Z, m_A;
  const Nuclear::Parameters m_nuc_params;
  std::unique_ptr<HF::HartreeFock> m_pHF{nullptr};

public:
  //! Sigma, correlation potential (2nd order)
  std::unique_ptr<MBPT::CorrelationPotential> m_Sigma{nullptr};

public:
  //! Nuclear potential
  std::vector<double> vnuc{};
  //! Direct/local part of the electron potential
  std::vector<double> vdir{};
  //! QED radiative potential
  RadiativePotential::Vrad vrad{};

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
  //! Number of electrons in the core
  int Ncore() const { return num_core_electrons; }
  const Nuclear::Parameters &get_nuclearParameters() const {
    return m_nuc_params;
  }
  bool exclude_exchangeQ() const {
    if (m_pHF == nullptr)
      return true;
    return m_pHF->excludeExchangeQ();
  }

  // Kill this (used once?)
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
  //! e.g., "Cs"
  std::string atomicSymbol() const { return AtomData::atomicSymbol(m_Z); }

  //! Prints table of core orbitals + energies etc. Optionally sorted by energy
  void printCore(bool sorted = true) const;

  //! @brief Prints table of valence orbitals + energies etc. Optionally sorted
  //! by energy
  //! @details Can optionally give it any list of orbitals to print
  void printValence(bool sorted = true,
                    const std::vector<DiracSpinor> &tmp_orbitals = {}) const;
  //! Prints table of Basis orbitals, compares to HF orbitals
  void printBasis(bool sorted = false) const;
  void printSpectrum(bool sorted = false) const;
  bool isInCore(int n, int k) const;
  bool isInValence(int n, int k) const;
  bool isInCore(const DiracSpinor &phi) const; // kill this one
  //! Largest n for core states
  int maxCore_n(int ka_in = 0) const;
  //! Largest l for core states
  int maxCore_l() const;

  //! Calculated rho(r) = sum_c psi^2(r) for core states
  std::vector<double> coreDensity() const;

  //! Performs hartree-Fock procedure for core: note: poplulates core
  void hartreeFockCore(HF::Method method, const std::string &in_core,
                       double eps_HF = 0, double h_d = 0, double g_t = 0);

  //! Calculates HF core energy (doesn't include magnetic QED?)
  auto coreEnergyHF() const;

  //! Performs hartree-Fock procedure for valence: note: poplulates valnece
  void hartreeFockValence(const std::string &in_valence_str,
                          const bool print = true);
  //! Forms Bruckner valence orbitals: (H_hf + Sigma)|nk> = e|nk>.
  void hartreeFockBrueckner(const bool print = true);
  //! First, fits Sigma to energies, then forms fitted Brueckner orbitals
  void fitSigma_hfBrueckner(const std::string &valence_list,
                            const std::vector<double> &fit_energies);
  //! Second-order MBPT energy shifts, calculates + prints
  void SOEnergyShift();

  //! Calculates radiative potential. Stores in vnuc, and Hmag
  void radiativePotential(double x_simple, double x_Ueh, double x_SEe_h,
                          double x_SEe_l, double x_SEm, double rcut,
                          double scale_rN, const std::vector<double> &x_spd);

  //! Calculates + populates basis [see BSplineBasis]
  void formBasis(const std::string &states_str, const std::size_t n_spl,
                 const std::size_t k_spl, const double r0_spl,
                 const double r0_eps, const double rmax_spl,
                 const bool positronQ = false);

  //! Calculates + populates Spectrum [see BSplineBasis]
  void formSpectrum(const std::string &states_str, const std::size_t n_spl,
                    const std::size_t k_spl, const double r0_spl,
                    const double r0_eps, const double rmax_spl,
                    const bool positronQ = false);

  //! Forms + stores correlation potential Sigma
  void formSigma(const int nmin_core = 1, const bool form_matrix = true,
                 const std::vector<double> &lambdas = {});

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

  //! @brief Returns [min,max] r values for which the core density (given l) is
  //! larger than cutoff (= eps*max_value)
  //! @details Returns the r values (au) for which the value of rho =
  //! \sum|psi^2|(r) drops below cutoff. Sum goes over all m for given l.
  //! Cut-off defined as eps*max, where max is maximum value for rho(r). Returns
  //! for each l in the core.
  std::tuple<double, double> lminmax_core_range(int l, double eps = 0.0) const;

  //! Local potential, e.g., Vl = Vnuc + Vdir + Vrad_el(l) - can be l-dependent
  std::vector<double> get_Vlocal(int l = 0) const;
  const std::vector<double> &get_Hmag(int l = 0) const;

private:
  void determineCore(std::string str_core_in);
  // not static, since skip_core. Make static version??
  std::vector<AtomData::DiracSEnken>
  listOfStates_nk(int num_val, int la, int lb = 0, bool skip_core = true) const;
  static std::vector<std::size_t>
  sortedEnergyList(const std::vector<DiracSpinor> &tmp_orbs,
                   bool do_sort = false);
};
