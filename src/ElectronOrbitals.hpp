#pragma once
#include "DiracSpinor.hpp"
#include "Grid.hpp"
#include <string>
#include <vector>

enum class NucleusType { Fermi, spherical, zero };

static bool dummy_bool{};

//******************************************************************************
class ElectronOrbitals {

public:
  ElectronOrbitals(int in_z, int in_a, int in_ngp, double rmin, double rmax,
                   double var_alpha = 1);

public:
  // orbitals:
  std::vector<DiracSpinor> core_orbitals;
  std::vector<DiracSpinor> valence_orbitals;
  // std::vector<DiracSpinor> basis;    // XXX break into core+valence?

  const Grid rgrid;

  // Potentials
  std::vector<double> vnuc;
  std::vector<double> vdir; // direct/local part of the electron potential

private:
  // store internal value for alpha (allows variation)
  const double m_alpha;
  // Atom info:
  const int m_Z, m_A;
  // nuclus info:
  double m_c, m_t;

  // number of electrons in each core shell (non-rel)
  std::vector<int> num_core_shell; // XXX This is dumb - try to fix!?
  int num_core_electrons = 0;      // Nc = N - M
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
  };

  std::size_t getStateIndex(int n, int k, bool &is_valence = dummy_bool) const;

  std::string coreConfiguration() const { return m_core_string; }
  std::string coreConfiguration_nice() const {
    return AtomInfo::niceCoreOutput(m_core_string);
  }
  std::string nuclearParams() const;
  std::string atom() const {
    return AtomInfo::atomicSymbol(m_Z) + ", Z=" + std::to_string(m_Z) +
           " A=" + std::to_string(m_A);
  }
  void printCore(bool sorted = true) const;
  void printValence(bool sorted = true,
                    const std::vector<DiracSpinor> &tmp_orbitals = {}) const;
  bool isInCore(int n, int k) const;
  int maxCore_n(int ka_in = 0) const;

  std::vector<std::size_t>
  sortedEnergyList(const std::vector<DiracSpinor> &tmp_orbs,
                   bool do_sort = false) const;

  // re-write this in terms of nkens !! XXX
  std::vector<std::vector<int>> listOfStates_nk(int num_val, int la, int lb = 0,
                                                bool skip_core = true) const;

public:
  void formNuclearPotential(NucleusType nucleus_type, double rc = 0,
                            double t = 0);

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
