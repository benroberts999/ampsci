#pragma once
#include "IO/FRW_fileReadWrite.hpp"
#include "LinAlg/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <array>
#include <iostream>
#include <optional>
#include <string_view>
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
  @brief jj -> LS recoupling amplitude for an antisymmetrised two-electron CSF.
  @details
  Returns the amplitude of the antisymmetrised jj-coupled CSF
  \f$ |\{(n_1 l_1 j_1)(n_2 l_2 j_2)\}; J\rangle \f$ (orbitals in stored,
  i.e., sorted, order) onto the antisymmetrised LS-coupled state
  \f$ |\{(n_1 l_1)(n_2 l_2)\} L S; J\rangle \f$ of the same non-relativistic
  configuration:

  \f[
    A(L,S) = \eta \sqrt{[j_1][j_2][L][S]}
    \begin{Bmatrix} l_1 & l_2 & L \\ 1/2 & 1/2 & S \\ j_1 & j_2 & J \end{Bmatrix}
  \f]

  Taken in the non-relativistic limit: the radial orbitals of
  \f$ j = l \pm 1/2 \f$ are treated as identical (overlap = 1).

  For a common non-relativistic shell (\f$ n_1 = n_2 \f$, \f$ l_1 = l_2 \f$)
  only L+S even terms exist (Pauli). When additionally \f$ j_1 \neq j_2 \f$
  the L+S odd components cancel in the antisymmetrisation and the even ones
  carry \f$ \eta = \sqrt{2} \f$; otherwise \f$ \eta = 1 \f$.
  In all cases \f$ \sum_{LS} A^2 = 1 \f$.

  @param n1,l1,twoj1  Quantum numbers of the first stored orbital.
  @param n2,l2,twoj2  Quantum numbers of the second stored orbital.
  @param L,S          Total orbital and spin angular momenta of the LS term.
  @param twoJ         Twice the total angular momentum 2J.
  @return Recoupling amplitude A(L,S); zero if forbidden.

  @note The sign convention follows the stored (sorted) orbital order; since
        nk_index sorting keeps the (n, l) order identical for all CSFs of one
        non-relativistic configuration, relative signs between such CSFs are
        consistent.
*/
double LS_amplitude(int n1, int l1, int twoj1, int n2, int l2, int twoj2, int L,
                    int S, int twoJ);

/*!
  @brief Expectation values of L^2 and S^2 for a two-electron CI state.
  @details
  Recouples each CSF to LS coupling (see @ref LS_amplitude) and accumulates,
  per non-relativistic configuration g,
  \f$ B_g(L,S) = \sum_{I \in g} c_I A_I(L,S) \f$, giving

  \f[
    \langle L^2 \rangle = \sum_{g,L,S} B_g(L,S)^2 \, L(L+1), \qquad
    \langle S^2 \rangle = \sum_{g,L,S} B_g(L,S)^2 \, S(S+1).
  \f]

  These are expectation values of the CI state, not eigenvalues: deviation
  from L(L+1), S(S+1) measures the LS-purity of the state.

  @param coefs  CI expansion coefficients (one per CSF).
  @param csfs   The CSF basis (matching @p coefs).
  @param twoJ   Twice the total angular momentum 2J.
  @return Pair {<L^2>, <S^2>}.

  @note Non-relativistic limit: radial overlaps between j = l +- 1/2 orbitals
        are set to 1, so for a normalised state the total LS weight is exactly
        1 and no renormalisation is required.
*/
std::pair<double, double>
expectation_L2S2(const LinAlg::View<const double> &coefs,
                 const std::vector<CSF2> &csfs, int twoJ);

//==============================================================================
/*!
  @brief Identifies one CI level: its (J, parity), and which solution.
  @details
  The standard text form is `J{+,-}:index`, e.g., `2+:3` is the fourth solution
  (index counts from zero) of the J=2, even parity, CI. 
  The index may be omitted, in which case it is zero:
  `0+` is the lowest even-parity J=0 solution. See @ref parse_level and
  @ref to_string.
*/
struct Level {
  //! Twice the total angular momentum, 2J
  int twoJ{0};
  //! Parity: +1 or -1
  int parity{1};
  //! Which solution, counting from zero, in order of energy
  std::size_t index{0};
};

/*!
  @brief Parses the text form of a CI level reference; see @ref Level.
  @details
  Accepts `2+:3` (standard) and `e2:3`; the index is optional. Surrounding
  whitespace is ignored. Only integer J is accepted, since only two-electron CI
  is implemented.

  @param str  Text form, e.g., `2+:3`, `e2:3`, `0+`.
  @return The level; empty if @p str is not a valid level reference.
*/
[[nodiscard]] std::optional<Level> parse_level(std::string_view str);

//! Text form of a CI level reference, e.g., "2+:3"; see @ref Level
[[nodiscard]] std::string to_string(const Level &level);

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
  //! Expectation value of L^2 for the CI state (-1 if not computed)
  double L2{-1.0};
  //! Expectation value of S^2 for the CI state (-1 if not computed)
  double S2{-1.0};
};

//==============================================================================
/*!
  @brief Container for CI solutions in a single J/Pi.
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

  /*!
    @brief Stores a single solution directly, without diagonalising.
    @details
    Replaces any existing solutions with the single one given.  For states that
    are not eigenstates of the CI Hamiltonian: e.g., the mixed states of
    @ref solve_mixed_state, for which @p energy is that of the reference state.

    @param energy Energy to be associated with the solution (atomic units).
    @param coefs  CI expansion coefficients; one per CSF.
  */
  void set_solution(double energy, const LinAlg::Vector<double> &coefs);

  //! Set configuration info for the ith solution (must be called manually after solve())
  void update_config_info(std::size_t i, const ConfigInfo &info);

  //! Full list of CSFs spanning this J/Pi
  const std::vector<CSF2> &CSFs() const;

  //! Returns reference to the ith CSF
  const CSF2 &CSF(std::size_t i) const;

  //! Energy of the ith CI solution (atomic units)
  double energy(std::size_t i) const;

  //! CI expansion coefficients for the ith solution (one per CSF)
  LinAlg::View<const double> coefs(std::size_t i) const;

  //! CI coefficient for the ith solution corresponding to the jth CSF
  double coef(std::size_t i, std::size_t j) const;

  //! Parity of the J/Pi (+/-1)
  int parity() const;

  //! Twice the total angular momentum 2J for this J/Pi
  int twoJ() const;

  //! Number of CI solutions currently stored
  std::size_t num_solutions() const;

  //! Configuration info for the ith solution (must have been set via update_config_info())
  const ConfigInfo &info(std::size_t i) const;

  /*!
    @brief Reads or writes CI solutions (energies, eigenvectors) to/from a multi-block binary file.
    @details
    A single file holds multiple blocks, one per (twoJ, parity) pair.
    Each block is self-describing: (twoJ, pi, num_csfs, num_solutions, E[num_csfs], M[num_csfs x num_csfs]).
    Energies and eigenvectors are always stored at full num_csfs size (zero-padded
    if only a partial solve was done), so block size is determined from num_csfs
    alone, enabling O(N) scan and in-place overwrite.

    The CSF basis (m_CSFs) is not touched; it must be constructed via the normal
    constructor before calling this function. The basis (and hence num_csfs) must
    be consistent with the file on write; a mismatch causes failure.

    ConfigInfo is not stored -- call update_config_info() after reading if needed.

    On read: scans blocks until matching (twoJ, pi) is found, verifies num_csfs,
    reads num_solutions, E, and M; returns false if block not found or num_csfs
    mismatches.

    On write: if the block already exists and num_csfs matches, overwrites it
    in-place using @ref IO::FRW::update. If it is a new block, appends it to
    the end of the file -- no rewrite of existing data.

    @note No settings are stored in the file: the filename identifies the
    calculation. The default filenames encode the settings that change the
    solutions -- see @ref CI::configuration_interaction.

    @param fname      Path to the binary file.
    @param rw         @ref IO::FRW::read to read; @ref IO::FRW::write to write.
    @param outstream  Output stream for progress/notes.

    @return True on success; false if the file does not exist (read), the block
            is not found (read), or num_csfs mismatches.
  */
  bool read_write(const std::string &fname, IO::FRW::RoW rw,
                  std::ostream &outstream = std::cout);
};

} // namespace CI
