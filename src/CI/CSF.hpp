#pragma once
#include "LinAlg/LinAlg.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <array>
#include <utility>
#include <vector>

namespace CI {

//==============================================================================
//! Very basic two-electron CSF. Only two-electron is implemented.
class CSF2 {
  int m_parity;

public:
  // nb: array of states is always sorted
  std::array<DiracSpinor::Index, 2> states;

  CSF2(const DiracSpinor &a, const DiracSpinor &b);

  DiracSpinor::Index state(std::size_t i) const;

  friend bool operator==(const CSF2 &A, const CSF2 &B);
  friend bool operator!=(const CSF2 &A, const CSF2 &B);

  //! Returns number of different orbitals between two CSFs
  static int num_different(const CSF2 &A, const CSF2 &B);

  //! returns _different_ orbitals, for case where CSFs differ by 1.
  //! i.e., returns {n,a} where |V> = |X_a^n> (i.e., V has n, but not a)
  static std::array<DiracSpinor::Index, 2> diff_1_na(const CSF2 &V,
                                                     const CSF2 &X);
  //! Returns the state in A and B that is the same (assumes A and B differ by 1)
  static DiracSpinor::Index same_1_j(const CSF2 &A, const CSF2 &B);

  //! Parity of the CSF, +/-1
  int parity() const;

  //! Single-particle configuration as a string, in relativistic or non-rel form
  std::string config(bool relativistic = false) const;
};

//==============================================================================
//! Forms list of all possible (2-particle) CSFs with given J and parity
std::vector<CSF2> form_CSFs(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis);

//==============================================================================
//! Basic configuration info for each CI level solution
struct ConfigInfo {
  // configuration (typically uses non-rel notation)
  std::string config{};
  // square of config coeficient (if non-rel, sum of all)
  double ci2{0.0};
  double gJ{0.0};
  double L{-1.0};
  double twoS{-1.0};
};

//==============================================================================
//! Stores the CI Solutions for given J and parity (only two-electron).
class PsiJPi {

  int m_twoj{-1};
  int m_pi{0};

  // Number of solutions stored:
  std::size_t m_num_solutions{0};
  // List of CSFs
  std::vector<CSF2> m_CSFs{};
  // Energy, and CI expansion coeficients
  std::pair<LinAlg::Vector<double>, LinAlg::Matrix<double>> m_Solution{};
  std::vector<ConfigInfo> m_Info{};

public:
  //! Construct containter for CI solutions.
  //! Constructs the CSFs, but doesn't solve system; have to call solve().
  PsiJPi(int twoJ, int pi, const std::vector<DiracSpinor> &cisp_basis)
      : m_twoj(twoJ), m_pi(pi), m_CSFs(form_CSFs(twoJ, pi, cisp_basis)) {}

  PsiJPi() {}

  //! Solves the CI equation for Hamiltonian matrix Hci;
  //! finds first num_solutions solutions.
  //! Doesn't set the Config info: have to call update_config_info() manually.
  //! @details If num_solutions=0 (or not given), will calculate _all_ solutions
  void solve(const LinAlg::Matrix<double> &Hci, int num_solutions = 0);

  //! You must manuall update the config. info for each solution (if required)
  void update_config_info(std::size_t i, const ConfigInfo &info);

  //! Full list of CSFs
  const std::vector<CSF2> &CSFs() const;

  //! The ith CSF
  const CSF2 &CSF(std::size_t i) const;

  //! Energy of the ith CI solution
  double energy(std::size_t i) const;

  //! List of CI expansion coefs for the ith CI solution
  LinAlg::View<const double> coefs(std::size_t i) const;

  //! The CI coeficient for the ith CI solution, corresponding to the jth CSF
  double coef(std::size_t i, std::size_t j) const;

  //! Parity for the CI solutions (+/-1)
  int parity() const;

  //! 2J for the CI solutions
  int twoJ() const;

  //! Number of CI solutions stored
  std::size_t num_solutions() const;

  //! Configuration info for the ith CI solution, if it has been set
  const ConfigInfo &info(std::size_t i) const;
};

} // namespace CI