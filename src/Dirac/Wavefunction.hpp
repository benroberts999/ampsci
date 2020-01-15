#pragma once
#include "Dirac/BSplineBasis.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "HF/HartreeFockClass.hpp" // forward decl..
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp" // NonRelSEConfig
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp" //PhysConst::alpha
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static bool dummy_bool{};
//******************************************************************************
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
  std::vector<DiracSpinor> core_orbitals;
  std::vector<DiracSpinor> valence_orbitals;
  std::vector<DiracSpinor> basis;
  const Grid rgrid;

private:
  const double m_alpha; // store internal value for alpha (allows variation)
  const int m_Z, m_A;
  Nuclear::Parameters m_nuc_params;
  std::unique_ptr<HartreeFock> m_pHF = nullptr;

public:
  // const
  std::vector<double> vnuc;
  std::vector<double> vdir;    // direct/local part of the electron potential
  std::vector<double> Hse_mag; // magnetic form-factor

private:
  // Core configuration (non-rel terms)
  std::vector<NonRelSEConfig> m_core_configs;
  int num_core_electrons = 0; // Nc = N - M
  std::string m_core_string = "";

public: // const methods: "views" into WF object
  // Rule is: if function is single-line, define here. Else, in .cpp
  double get_alpha() const { return m_alpha; }
  int Znuc() const { return m_Z; }
  int Anuc() const { return m_A; }
  int Nnuc() const { return (m_A > m_Z) ? (m_A - m_Z) : 0; }
  int Ncore() const { return num_core_electrons; }
  const Nuclear::Parameters &get_nuclearParameters() const {
    return m_nuc_params;
  }
  bool exclude_exchangeQ() const {
    if (m_pHF == nullptr)
      return true;
    return m_pHF->m_excludeExchange;
  }

  auto get_VexPsi(const DiracSpinor &psi) const {
    // XXX add check!? XXX
    return m_pHF->vex_psia(psi);
  }

  std::size_t getStateIndex(int n, int k, bool &is_valence = dummy_bool) const;
  const DiracSpinor *getState(int n, int k,
                              bool &is_valence = dummy_bool) const;

  std::string coreConfiguration() const { return m_core_string; }
  std::string coreConfiguration_nice() const {
    return AtomData::niceCoreOutput(m_core_string);
  }
  std::string nuclearParams() const;

  std::string atom() const {
    return AtomData::atomicSymbol(m_Z) + ", Z=" + std::to_string(m_Z) +
           " A=" + std::to_string(m_A);
  }
  std::string atomicSymbol() const { return AtomData::atomicSymbol(m_Z); }
  void printCore(bool sorted = true) const;
  void printValence(bool sorted = true,
                    const std::vector<DiracSpinor> &tmp_orbitals = {}) const;
  void printBasis(bool sorted = false) const;
  bool isInCore(int n, int k) const;
  bool isInValence(int n, int k) const;
  bool isInCore(const DiracSpinor &phi) const;
  int maxCore_n(int ka_in = 0) const;
  int maxCore_l() const;

  // not static, since skip_core. Make static version??
  std::vector<DiracSEnken> listOfStates_nk(int num_val, int la, int lb = 0,
                                           bool skip_core = true) const;

  std::vector<double> coreDensity() const;

  static std::vector<std::size_t>
  sortedEnergyList(const std::vector<DiracSpinor> &tmp_orbs,
                   bool do_sort = false);

public:
  void solveDirac(DiracSpinor &psi, double e_a, const std::vector<double> &vex,
                  int log_dele_or = 0) const;
  void solveDirac(DiracSpinor &psi, double e_a = 0, int log_dele_or = 0) const;

  void solveInitialCore(const std::string &str_core_in, int log_dele_or = 0);
  void solveNewValence(int n, int k, double en_a = 0, int log_dele_or = 0);

  static void orthonormaliseOrbitals(std::vector<DiracSpinor> &tmp_orbs,
                                     int num_its = 1);
  static void orthonormaliseWrt(DiracSpinor &psi_v,
                                const std::vector<DiracSpinor> &in_orbs);
  static void orthogonaliseWrt(DiracSpinor &psi_v,
                               const std::vector<DiracSpinor> &in_orbs);

  void hartreeFockCore(HFMethod method, const std::string &in_core,
                       double eps_HF = 0, double h_d = 0, double g_t = 0);
  auto coreEnergyHF() const;
  void hartreeFockValence(const std::string &in_valence_str);

  void radiativePotential(double x_Ueh, double x_SEe_h, double x_SEe_l,
                          double x_SEm, double rcut, double scale_rN);

  double enGuessCore(int n, int l) const;
  double enGuessVal(int n, int ka) const;

  void formBasis(const std::string &states_str, const std::size_t n_spl,
                 const std::size_t k_spl, const double r0_spl,
                 const double rmax_spl);

private:
  void determineCore(std::string str_core_in);
};
