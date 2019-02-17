#pragma once
#include "DiracSpinor.h"
#include "Grid.h"
#include <string>
#include <vector>

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
  // std::vector<std::vector<double>> f;
  // std::vector<std::vector<double>> g;
  // std::vector<double> en;
  // state info:
  // std::vector<int> nlist;
  // std::vector<int> kappa;
  // info from solveing DE
  // std::vector<int> pinflist;   // practical infinity
  // std::vector<int> itslist;    // num. iterations
  // std::vector<double> epslist; // convergance
  // occupancy fraction.
  // Note: avg over non-rel for core, but rel for valence!
  // std::vector<double> occ_frac;

  std::vector<int> stateIndexList;
  std::vector<int> coreIndexList;
  std::vector<int> valenceIndexList;

  int num_core_states = 0;
  int num_core_electrons = 0; // Nc = N - M

  const Grid rgrid;

  // Potentials
  std::vector<double> vnuc;
  std::vector<double> vdir; // direct/local part of the electron potential

private:
  // store internal value for alpha (allows variation)
  const double alpha;
  // Atom info:
  const int Z_, A_;
  // number of electrons in each core shell (non-rel??)
  std::vector<int> num_core_shell;

public:
  double get_alpha() const;
  int Znuc() const;
  int Anuc() const;
  int Nnuc() const;
  int ka(int i) const;
  int n_pqn(int i) const;
  int lorb(int i) const;
  int twoj(int i) const;
  double jtot(int i) const;
  double rinf(int i) const;
  int getStateIndex(int n, int k, bool forceVal = false) const;
  int getRadialIndex(double r_target) const;
  int numberOfStates() const; // statesCount() ?

  std::string seTermSymbol(int ink, bool gnuplot = false) const;

  double radialIntegral(int a, int b, Operator op = Operator::unity) const;
  double radialIntegral(int a, int b, const std::vector<double> &vint,
                        Operator op = Operator::unity) const;

  int solveLocalDirac(int n, int k, double en_a = 0, int log_dele_or = 0,
                      bool iscore = false);
  int reSolveDirac(unsigned long i, double e_a = 0, int log_dele_or = 0);
  int reSolveDirac(unsigned long i, double e_a, const std::vector<double> &vex,
                   int log_dele_or = 0);

  void formNuclearPotential(NucleusType nucleus_type, double rc = 0,
                            double t = 0);

  double diracen(double z, double n, int k) const;

  int solveInitialCore(std::string str_core_in, int log_dele_or = 0);
  bool isInCore(int n, int k) const;
  int maxCore_n(int ka_in = 0) const;

  void orthonormaliseOrbitals(int num_its = 1);
  void orthonormaliseValence(int iv, int num_its);

  void sortedEnergyList(std::vector<int> &sort_list,
                        bool do_sort = false) const;

private:
  int determineCore(std::string str_core_in);

  void sphericalNucleus(double rnuc = 0);
  void fermiNucleus(double t = 0, double c = 0);
  void zeroNucleus();

  double enGuessCore(int n, int l) const;
  double enGuessVal(int n, int ka) const;
};
