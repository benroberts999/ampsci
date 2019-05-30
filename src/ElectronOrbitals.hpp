#pragma once
#include "DiracSpinor.hpp"
#include "Grid.hpp"
#include <string>
#include <vector>
// class DiracSpinor;
// class Grid;

enum class Operator { unity, gamma0, gamma5, dr, dr2 };
enum class NucleusType { Fermi, spherical, zero };

//******************************************************************************
class ElectronOrbitals {

public:
  ElectronOrbitals(int in_z, int in_a, int in_ngp, double rmin, double rmax,
                   double var_alpha = 1);

public:
  // orbitals:
  std::vector<DiracSpinor> orbitals;

  const Grid rgrid;

  // Potentials
  std::vector<double> vnuc;
  std::vector<double> vdir; // direct/local part of the electron potential

  std::vector<std::size_t> stateIndexList;
  std::vector<std::size_t> coreIndexList;
  std::vector<std::size_t> valenceIndexList;

private:
  // store internal value for alpha (allows variation)
  const double m_alpha;
  // Atom info:
  const int m_Z, m_A;

  // number of electrons in each core shell (non-rel??)
  std::vector<int> num_core_shell;
  int num_core_electrons = 0; // Nc = N - M
  std::string m_core_string = "";

public:
  double get_alpha() const { return m_alpha; }
  int Znuc() const { return m_Z; }
  int Anuc() const { return m_A; }
  int Nnuc() const { return (m_A > m_Z) ? (m_A - m_Z) : 0; }
  int Ncore() const { return num_core_electrons; }
  double rinf(const DiracSpinor &phi) const;
  int getRadialIndex(double r_target) const;
  std::size_t getStateIndex(int n, int k) const;

  std::string coreConfiguration() { return m_core_string; }
  std::string coreConfiguration_nice() {
    return ATI::niceCoreOutput(m_core_string);
  }

  // // XXX Combine with default argument?
  // double radialIntegral(const DiracSpinor &psi_a, const DiracSpinor &psi_b,
  //                       const std::vector<double> &vint,
  //                       Operator op = Operator::unity) const;
  // double radialIntegral(const DiracSpinor &psi_a, const DiracSpinor &psi_b,
  //                       Operator op = Operator::unity) const;

  // XXX ? best way?
  // XXX Should really have just 1! ?
  // Then, different way to contruct new spinor??
  void solveLocalDirac(int n, int k, double en_a = 0, int log_dele_or = 0,
                       bool iscore = false);
  void reSolveDirac(DiracSpinor &psi_a, double e_a = 0, int log_dele_or = 0);
  void reSolveDirac(DiracSpinor &psi_a, double e_a,
                    const std::vector<double> &vex, int log_dele_or = 0);

  void formNuclearPotential(NucleusType nucleus_type, double rc = 0,
                            double t = 0);

  int solveInitialCore(std::string str_core_in, int log_dele_or = 0);
  bool isInCore(int n, int k) const;
  int maxCore_n(int ka_in = 0) const;

  void orthonormaliseOrbitals(int num_its = 1);
  void orthonormaliseValence(DiracSpinor &psi_v, int num_its = 1,
                             bool core_only = false) const; // "const" OK?

  std::vector<int> sortedEnergyList(bool do_sort = false) const;

private:
  void determineCore(const std::string &str_core_in);

  double enGuessCore(int n, int l) const;
  double enGuessVal(int n, int ka) const;
};
