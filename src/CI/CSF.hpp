#pragma once
#include "LinAlg/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <array>
#include <optional>
#include <utility>
#include <vector>

namespace CI {

//==============================================================================
/*!
  @brief Two-electron configuration state function (CSF).
  @details
  A CSF is an antisymmetrised two-electron basis state with definite total
  angular momentum (\f$ J^2 \f$, \f$ J_z \f$) and parity, 
  built from a pair of single-particle
  relativistic orbitals.  Only two-electron CSFs are implemented.

  Each CSF2 stores the indices of its two constituent orbitals (always sorted
  to avoid double-counting) and the total parity, which is the product of the
  parities of the two single-particle states.

  @note The orbital pair is stored as a sorted array of DiracSpinor::Index
        (uint16_t) rather than DiracSpinor references, so CSF2 objects are
        cheap to copy and store. There is a limit to maximum n<=256 - see @ref Angular::nk_to_index
*/
class CSF2 {
  int m_parity;

public:
  // nb: array of states is always sorted
  std::array<DiracSpinor::Index, 2> states;

  CSF2(const DiracSpinor &a, const DiracSpinor &b);

  //! Index (nk_index) of the ith constituent orbital (i = 0 or 1)
  DiracSpinor::Index state(std::size_t i) const;

  friend bool operator==(const CSF2 &A, const CSF2 &B);
  friend bool operator!=(const CSF2 &A, const CSF2 &B);

  /*!
    @brief Returns the number of orbitals that differ between two CSFs (0, 1,
    or 2).
    @details
    Used to select the appropriate Slater-Condon rule when evaluating CI matrix
    elements: 0 -- diagonal; 1 -- single substitution; 2 -- double
    substitution; >2 -- zero by orthogonality.
  */
  static int num_different(const CSF2 &A, const CSF2 &B);

  /*!
    @brief For two CSFs differing by exactly one orbital, returns {n, a} where
    @p V contains orbital n and @p X contains orbital a.
    @details
    Identifies the "particle" index n (in @p V but not @p X) and the "hole"
    index a (in @p X but not @p V), as needed to apply the single-substitution
    Slater-Condon rule: \f$ \langle V | \hat{O} | X \rangle \f$ where
    \f$ |V\rangle = \hat{a}^\dag_n \hat{a}_a |X\rangle \f$.

    @warning Result is undefined if @p V and @p X do not differ by exactly one
             orbital; check with num_different() first.
  */
  static std::array<DiracSpinor::Index, 2> diff_1_na(const CSF2 &V,
                                                     const CSF2 &X);

  /*!
    @brief Returns the orbital index shared by two CSFs that differ by exactly
    one orbital.
    @details
    Extracts the common (spectator) orbital needed for single-substitution
    matrix elements.

    @warning Assumes @p A and @p B differ by exactly one orbital.
  */
  static DiracSpinor::Index same_1_j(const CSF2 &A, const CSF2 &B);

  //! Parity of the CSF, +/-1
  int parity() const;

  //! Single-particle configuration as a string, in relativistic or non-rel form
  std::string config(bool relativistic = false) const;
};

//==============================================================================
/*!
  @brief Forms all two-electron CSFs with given total J and parity.
  @details
  Iterates over all pairs of single-particle states in @p cisp_basis and
  retains those whose angular momenta can be coupled to total \f$ J = \f$
  @p twoJ /2 and whose combined parity equals @p parity.  Duplicate pairs are
  excluded by construction.

  @param twoJ       Twice the total angular momentum 2J.
  @param parity     Total parity: +1 (even) or -1 (odd).
  @param cisp_basis Single-particle basis from which CSFs are constructed.
  @return Sorted list of all valid two-electron CSFs for the given J and parity.
*/
std::vector<CSF2> form_CSFs(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis);

//==============================================================================
/*!
  @brief Configuration metadata for a single CI level.
  @details
  Stores identifying information derived after solving the CI eigenvalue
  problem: the dominant non-relativistic configuration label, the squared CI
  coefficient of that configuration, and approximate good quantum numbers
  (g_J factor, L, S) where they can be assigned.

  Fields are left at their default (empty/negative) values if not yet computed;
  call @ref PsiJPi::update_config_info() to populate them.
*/
struct ConfigInfo {
  //! Dominant configuration label (typically non-relativistic notation)
  std::string config{};
  //! Squared CI coefficient of the dominant configuration (or sum over non-rel degenerates)
  double ci2{0.0};
  double gJ{0.0};
  //! Approximate orbital angular momentum L (-1 if not assigned)
  double L{-1.0};
  //! Twice the approximate spin S (-1 if not assigned)
  double twoS{-1.0};
};

//==============================================================================
/*!
  @brief Container for CI solutions in a single (J, parity) sector.
  @details
  Holds the complete set of configuration state functions and the results of
  the CI diagonalisation for a fixed total angular momentum J and parity.

  Construction builds the CSF basis via form_CSFs() but does not solve the
  eigenvalue problem; call solve() separately after constructing the CI
  Hamiltonian matrix.  Configuration labels are not set automatically -- call
  update_config_info() for each solution after solving.

  @note Only two-electron (two-particle) systems are supported.
*/
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
  /*!
    @brief Constructs the CSF basis for the given J and parity; does not solve.
    @details
    Calls form_CSFs() to build the list of two-electron CSFs.  The eigenvalue
    problem is not solved until solve() is called with the CI Hamiltonian.

    @param twoJ        Twice the total angular momentum 2J.
    @param pi          Total parity: +1 or -1.
    @param cisp_basis  Single-particle basis used to construct the CSFs.
  */
  PsiJPi(int twoJ, int pi, const std::vector<DiracSpinor> &cisp_basis)
    : m_twoj(twoJ), m_pi(pi), m_CSFs(form_CSFs(twoJ, pi, cisp_basis)) {}

  PsiJPi() {}

  /*!
    @brief Solves the CI eigenvalue problem for the given Hamiltonian matrix.
    @details
    Diagonalises @p Hci and stores the resulting eigenvalues and eigenvectors.
    Does not populate ConfigInfo; call update_config_info() separately.

    - If @p num_solutions > 0, finds only the lowest @p num_solutions eigenpairs.
    - If @p all_below is set, finds all eigenpairs with energy below that value
      (in cm^-1); @p num_solutions is then ignored.
    - If both are unset (or @p num_solutions <= 0), all eigenpairs are computed.

    @param Hci           CI Hamiltonian matrix in the CSF basis.
    @param num_solutions Number of lowest solutions to find [0 = all].
    @param all_below     If set, find all solutions below this energy (cm^-1).
  */
  void solve(const LinAlg::Matrix<double> &Hci, int num_solutions = 0,
             std::optional<double> all_below = {});

  //! Set configuration info for the ith solution (must be called manually after solve())
  void update_config_info(std::size_t i, const ConfigInfo &info);

  //! Full list of CSFs spanning this (J, parity) sector
  const std::vector<CSF2> &CSFs() const;

  //! Returns reference to the ith CSF
  const CSF2 &CSF(std::size_t i) const;

  //! Energy of the ith CI solution (atomic units)
  double energy(std::size_t i) const;

  //! CI expansion coefficients for the ith solution (one per CSF)
  LinAlg::View<const double> coefs(std::size_t i) const;

  //! CI coefficient for the ith solution corresponding to the jth CSF
  double coef(std::size_t i, std::size_t j) const;

  //! Parity of the sector (+/-1)
  int parity() const;

  //! Twice the total angular momentum 2J for this sector
  int twoJ() const;

  //! Number of CI solutions currently stored
  std::size_t num_solutions() const;

  //! Configuration info for the ith solution (must have been set via update_config_info())
  const ConfigInfo &info(std::size_t i) const;
};

} // namespace CI
