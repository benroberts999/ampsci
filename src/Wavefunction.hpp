#pragma once
#include "DiracSpinor.hpp"
#include "Grid.hpp"
#include "HartreeFockClass.hpp" // forward decl..
#include "Physics/Nuclear.hpp"
#include <memory>
#include <string>
#include <vector>

static bool dummy_bool{};

//******************************************************************************
class Wavefunction {

public:
  Wavefunction(int in_z, const GridParameters &gridparams,
               const Nuclear::Parameters &nuc_params, double var_alpha = 1);

public:
  // orbitals:
  std::vector<DiracSpinor> core_orbitals;
  std::vector<DiracSpinor> valence_orbitals;
  const Grid rgrid;

private:
  // store internal value for alpha (allows variation)
  const double m_alpha;
  // Atom info:
  const int m_Z, m_A; /*don't need A twice (its inside nucl params!)*/
  Nuclear::Parameters m_nuc_params;
  std::unique_ptr<HartreeFock> m_pHF = nullptr;

public:
  const std::vector<double> vnuc;
  std::vector<double> vdir; // direct/local part of the electron potential

private:
  // Core configuration (non-rel terms)
  std::vector<NonRelSEConfig> m_core_configs;
  int num_core_electrons = 0; // Nc = N - M
  std::string m_core_string = "";

public:
  // Rule is: if function is single-line, define here. Else, in .cpp
  double get_alpha() const { return m_alpha; }
  int Znuc() const { return m_Z; }
  int Anuc() const { return m_A; }
  int Nnuc() const { return (m_A > m_Z) ? (m_A - m_Z) : 0; }
  int Ncore() const { return num_core_electrons; }
  double rinf(const DiracSpinor &phi) const { return rgrid.r[phi.pinf]; };
  int getRadialIndex(double r_target) const {
    return (int)rgrid.getIndex(r_target, true);
  }

  void hartreeFockCore(HFMethod method, const std::string &in_core,
                       double eps_HF = 0, double h_d = 0, double g_t = 0) {
    // XXX Update this (and HF) so that it doesn't re-Create m_pHF
    // AND, so that can re-run core!
    m_pHF = std::make_unique<HartreeFock>(
        HartreeFock(method, *this, in_core, eps_HF, h_d, g_t));
  }
  auto coreEnergyHF() const {
    if (m_pHF == nullptr) {
      std::cerr
          << "WARNING 62: Cant call coreEnergyHF before hartreeFockCore\n";
      return 0.0;
    }
    return m_pHF->calculateCoreEnergy();
  }
  void hartreeFockValence(const std::string &in_valence_str) {
    if (m_pHF == nullptr) {
      std::cerr << "WARNING 62: Cant call hartreeFockValence before "
                   "hartreeFockCore\n";
      return;
    }
    auto val_lst = AtomInfo::listOfStates_nk(in_valence_str);
    for (const auto &nk : val_lst) {
      if (!isInCore(nk.n, nk.k))
        m_pHF->solveNewValence(nk.n, nk.k);
    }
  }
  auto get_VexPsi(const DiracSpinor &psi) const { // XXX add check!?
    return m_pHF->vex_psia(psi);
  }

  std::size_t getStateIndex(int n, int k, bool &is_valence = dummy_bool) const;
  std::size_t getStateIndex(const DiracSpinor &psi,
                            bool &is_valence = dummy_bool) const;

  std::string coreConfiguration() const { return m_core_string; }
  std::string coreConfiguration_nice() const {
    return AtomInfo::niceCoreOutput(m_core_string);
  }
  std::string nuclearParams() const;

  std::string atom() const {
    return AtomInfo::atomicSymbol(m_Z) + ", Z=" + std::to_string(m_Z) +
           " A=" + std::to_string(m_A);
  }
  std::string atomicSymbol() const { return AtomInfo::atomicSymbol(m_Z); }
  void printCore(bool sorted = true) const;
  void printValence(bool sorted = true,
                    const std::vector<DiracSpinor> &tmp_orbitals = {}) const;
  bool isInCore(int n, int k) const;
  bool isInCore(const DiracSpinor &phi) const;
  int maxCore_n(int ka_in = 0) const;
  int maxCore_l() const;

  std::vector<std::size_t>
  sortedEnergyList(const std::vector<DiracSpinor> &tmp_orbs,
                   bool do_sort = false) const;

  std::vector<DiracSEnken> listOfStates_nk(int num_val, int la, int lb = 0,
                                           bool skip_core = true) const;

public:
  void solveDirac(DiracSpinor &psi, double e_a, const std::vector<double> &vex,
                  int log_dele_or = 0) const;
  void solveDirac(DiracSpinor &psi, double e_a = 0, int log_dele_or = 0) const;

  void solveInitialCore(std::string str_core_in, int log_dele_or = 0);
  void solveNewValence(int n, int k, double en_a = 0, int log_dele_or = 0);

  static void orthonormaliseOrbitals(std::vector<DiracSpinor> &tmp_orbs,
                                     int num_its = 1);
  void orthonormaliseWrtCore(DiracSpinor &psi_v) const;

private:
  void determineCore(const std::string &str_core_in);

  double enGuessCore(int n, int l) const;
  double enGuessVal(int n, int ka) const;
};
