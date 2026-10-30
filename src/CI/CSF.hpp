#pragma once
#include "ConfigurationSpace.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "LinAlg/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <array>
#include <functional>
#include <iostream>
#include <optional>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace CI {

//------------------------------------------------------------------------------
/*!
  @brief One-electron configuration state function (CSF): a single orbital.

  @details
  The simplest CSF: one single-particle relativistic orbital, with
  \f$ J = j \f$ and the parity of the orbital. All CSFs of one (J, parity)
  therefore have the same \f$ \kappa \f$, and the CI Hamiltonian between them
  is just the one-body \f$ h_1 \f$ (there is no two-body part). Allows
  single-valence systems to be treated with the CI machinery.
  See @ref CSF2 and @ref CSF3 for two and three electrons, and @ref CSF.

  @note The orbital is stored as a one-element array of DiracSpinor::Index,
        as for @ref CSF2 and @ref CSF3.
*/
class CSF1 {
public:
  std::array<DiracSpinor::Index, 1> states;

private:
  int m_parity;

public:
  explicit CSF1(const DiracSpinor &a);

  //! Index (nk_index) of the orbital (i must be 0)
  DiracSpinor::Index state(std::size_t i) const;

  friend bool operator==(const CSF1 &A, const CSF1 &B);
  friend bool operator!=(const CSF1 &A, const CSF1 &B);

  //! Number of orbitals that differ between two CSFs: 0 or 1
  static int num_different(const CSF1 &A, const CSF1 &B);

  //! Parity of the CSF, +/-1
  int parity() const { return m_parity; }

  //! Number of electrons in the CSF
  std::size_t num_electrons() const { return 1; }
  //! Number of particles (excited orbitals) in the CSF
  std::size_t num_particles() const { return 1; }
  //! Number of holes (core orbitals removed) in the CSF
  std::size_t num_holes() const { return 0; }

  //! Single-particle configuration as a string, in relativistic or non-rel form
  std::string config(bool relativistic = false) const;
};

//------------------------------------------------------------------------------
/*!
  @brief Two-electron configuration state function (CSF).

  @details
  A CSF is an antisymmetrised two-electron basis state with definite total
  angular momentum (\f$ J^2 \f$, \f$ J_z \f$) and parity,
  built from a pair of single-particle
  relativistic orbitals. See @ref CSF1 for one electron, and @ref CSF3 for
  three.

  For orbitals v and w, the determinant coefficients (see @ref CSF) are
  simply Clebcsh-Gordan coefficients:
  \f[
    |vw; J J_z\rangle = \eta_{vw} \sum_{m_v m_w}
    C^{J J_z}_{j_v m_v, j_w m_w} \, |n_v\kappa_v m_v, n_w\kappa_w m_w\rangle,
  \f]
  where \f$ \eta_{vw} = 1/\sqrt{2} \f$ if v = w (both electrons in the same
  orbital), and 1 otherwise. 
  
  The orbitals are put in sorted order (v before w); this fixes the phase
  convention of all CSF2 matrix elements (e.g., CSF2_Coulomb()).

  Each CSF2 stores the indices (@ref DiracSpinor::Index) of its two constituent orbitals (always sorted to avoid double-counting) and the total parity, 
  which is the product of the parities of the two single-particle states.

  @note The orbital pair is stored as a sorted array of DiracSpinor::Index
        (uint16_t) rather than DiracSpinor references, so CSF2 objects are
        cheap to copy and store. There is a limit to maximum n<=256 - see @ref Angular::nk_to_index
*/
class CSF2 {

public:
  // nb: array of states is always sorted
  std::array<DiracSpinor::Index, 2> states;

private:
  int m_parity;

public:
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
  int parity() const { return m_parity; }

  //! Number of electrons in the CSF
  std::size_t num_electrons() const { return 2; }
  //! Number of particles (excited orbitals) in the CSF
  std::size_t num_particles() const { return 2; }
  //! Number of holes (core orbitals removed) in the CSF
  std::size_t num_holes() const { return 0; }

  //! Single-particle configuration as a string, in relativistic or non-rel form
  std::string config(bool relativistic = false) const;
};

//------------------------------------------------------------------------------
/*!
  @brief Three-electron configuration state function (CSF).

  @details
  An antisymmetrised three-electron state with definite total angular momentum
  J, J_z, and parity, built from three single-particle orbitals (which may
  coincide).

  The determinant coefficients (see @ref CSF) are not stored directly. Instead,
  the CSF is stored as an expansion over "parent" states, which are
  antisymmetric in electrons 1 and 2 only: a two-electron CSF \f$ [xy] \f$
  (exactly as @ref CSF2) with angular momentum \f$ J_2 \f$, coupled to the
  orbital z of the third electron:
  \f[
    |[xy] J_2, z; J J_z\rangle = \sum_{M_2 m_z}
    C^{J J_z}_{J_2 M_2, j_z m_z} \, |xy; J_2 M_2\rangle_{12} \,
    |n_z\kappa_z m_z\rangle_3 ,
  \f]

  where \f$ |xy; J_2 M_2\rangle_{12} \f$ is the CSF2 state of 
  electrons 1 and 2. 
  
  The three-electron CSF is then:
  
  \f[
    |\Psi; J J_z\rangle = \sum_p C_p \, |[x_p y_p] J_{2p}, z_p; J J_z\rangle,
  \f]

  normalised such that

  \f[
    \sum_p C_p^2 = 1.
  \f]

  The parent coefficients \f$ C_p \f$ are chosen such that
  \f$ |\Psi\rangle \f$ is fully antisymmetric (see @ref form_CSF3s); together
  with the Clebsch-Gordan coefficients above, they fix the determinant
  coefficients \f$ d_{\{m\}} \f$. (In the literature, these are called
  "coefficients of fractional parentage": the three-electron state is built
  by adding one electron to each two-electron parent, and \f$ C_p \f$ is the
  amplitude of parent p.)

  Expanding over parents reduces three-electron matrix elements to sums of
  two-electron (or one-electron) ones. E.g., for a two-body operator,
  \f[
    \langle\Psi|\sum_{i<j} g_{ij}|\Psi'\rangle
    = 3\sum_{pq} C_p C'_q \, \delta_{z_p z_q}\delta_{J_{2p}J_{2q}}
    \langle [x_p y_p] J_{2p}|g|[x_q y_q] J_{2q}\rangle .
  \f]

  The expansion depends on J, so (as for @ref CSF2) each CSF belongs to a
  single J. Formed by @ref form_CSF3s.

  @note Orbitals are stored as DiracSpinor::Index (sorted), as in @ref CSF2.
*/
class CSF3 {
public:
  //! One parent state of a CSF3, and its coefficient
  struct Parent {
    //! Orbitals of the parent two-electron CSF (sorted, as in CSF2)
    std::array<DiracSpinor::Index, 2> pair;
    //! Orbital of the third electron
    DiracSpinor::Index spectator;
    //! Twice the angular momentum of the parent pair, 2*J2
    int twoJ2;
    //! Coefficient of this parent, C_p
    double coef;
  };

  //! The three orbitals; always sorted
  std::array<DiracSpinor::Index, 3> states;

private:
  int m_twoJ12;
  std::vector<Parent> m_parents;
  int m_parity;

public:
  /*!
    @brief Constructs from the three orbitals and the parent expansion.
    @param orbitals  The three orbitals (any order; stored sorted).
    @param twoJ12    Label: 2*J12 of the coupling the CSF was formed from (see
                     @ref form_CSF3s).
    @param parents   Parent states and their coefficients.
  */
  CSF3(const std::array<DiracSpinor::Index, 3> &orbitals, int twoJ12,
       std::vector<Parent> parents);

  //! Index (nk_index) of the ith constituent orbital (i = 0, 1, or 2)
  DiracSpinor::Index state(std::size_t i) const;

  //! Label: 2*J12 of the two-electron coupling the CSF was formed from
  int twoJ12() const;

  //! Parent states and their coefficients
  const std::vector<Parent> &parents() const;

  /*!
    @brief Returns the number of orbitals that differ between two CSFs (0 to
    3), counting repeated orbitals.
    @details
    Matrix elements of one- and two-body operators vanish if more than two
    orbitals differ.
  */
  static int num_different(const CSF3 &A, const CSF3 &B);

  //! Parity of the CSF, +/-1
  int parity() const { return m_parity; }

  //! Number of electrons in the CSF
  std::size_t num_electrons() const { return 3; }
  //! Number of particles (excited orbitals) in the CSF
  std::size_t num_particles() const { return 3; }
  //! Number of holes (core orbitals removed) in the CSF
  std::size_t num_holes() const { return 0; }

  //! Single-particle configuration as a string, in relativistic or non-rel form
  std::string config(bool relativistic = false) const;

  friend bool operator==(const CSF3 &A, const CSF3 &B);
  friend bool operator!=(const CSF3 &A, const CSF3 &B);
};

//------------------------------------------------------------------------------
/*!
  @brief Exchange of electrons 2 and 3, between two parent states of a
  CSF3.
  @details
  Matrix element
  \f$ \langle [x'y']J_2', z'; J | P_{23} | [xy]J_2, z; J \rangle \f$, where
  each state is antisymmetric in electrons 1 and 2 (see @ref CSF3::Parent). In
  terms of uncoupled-pair states, the exchange is the recoupling

  \f[
    P_{23}|(o_1 o_2)J_{12}, o_3; J\rangle
    = \sum_K (-1)^{j_3+j_2+K+J_{12}}\sqrt{[K][J_{12}]}
      \begin{Bmatrix} j_3 & j_1 & K \\ j_2 & J & J_{12} \end{Bmatrix}
      |(o_1 o_3)K, o_2; J\rangle .
  \f]

  @param bra,ket  Parent states (their coef members are not used).
  @param twoJ     Twice the total angular momentum 2J.
  @return The matrix element.
*/
double parent_exchange_23(const CSF3::Parent &bra, const CSF3::Parent &ket,
                          int twoJ);

/*!
  @brief Forms all three-electron CSFs of a single configuration, with given
  total J.
  @details
  Builds an orthonormal basis for the antisymmetric three-electron states of
  the configuration, each as its expansion over parent states (see
  @ref CSF3):

  1. Parent basis: all \f$ |[xy]J_2, z; J\rangle \f$, for each distinct choice
     of orbital z for the third electron, with [xy] the remaining pair. J2 is
     even if x = y, and is limited by the triangle rules. These are
     orthonormal, and span the states antisymmetric in electrons 1 and 2.
  2. On that space, the full antisymmetriser is
     \f[
       A = \frac{1}{3}(1 - 2P_{23})
     \f]
     (see @ref parent_exchange_23); each column \f$ A|p\rangle \f$ is fully
     antisymmetric.
  3. Gram-Schmidt orthonormalisation of the \f$ A|p\rangle \f$, starting with
     the standard couplings ((ab)J12, c)J (or ((bc)J12, a)J if only b = c);
     linearly dependent vectors are discarded.

  For three distinct orbitals this gives one CSF per allowed J12; for
  \f$ j^2 j' \f$, one per even J12; for \f$ j^3 \f$, the standard
  \f$ j^3 \f$ states, including repeated J (e.g., the two J = 9/2 states of
  \f$ (9/2)^3 \f$).

  @param orbitals  The three orbitals (nk_index; any order; may coincide).
  @param twoJ      Twice the total angular momentum 2J.
  @return The CSFs of this configuration and J; empty if there are none.
*/
std::vector<CSF3> form_CSF3s(std::array<DiracSpinor::Index, 3> orbitals,
                             int twoJ);

//------------------------------------------------------------------------------
/*!
  @brief Two-particle, one-hole configuration state function (CSF).

  @details
  A hole excitation of the closed-shell core: two electrons in excited
  orbitals v and w, and one electron removed from the core orbital a. With
  the hole creation operator \f$ h^\dagger_{am} = (-1)^{j_a - m} a_{a,-m} \f$
  (a rank-\f$ j_a \f$ tensor), the CSF is the pair of particles (an
  antisymmetrised two-electron state, exactly as @ref CSF2, with angular
  momentum \f$ J_{12} \f$) coupled to the hole:
  \f[
    |(vw)J_{12}, \tilde a; J J_z\rangle = \eta_{vw} \sum
    C^{J_{12} M_{12}}_{j_v m_v, j_w m_w} C^{J J_z}_{J_{12} M_{12}, j_a m}
    \, a^\dagger_{v m_v} a^\dagger_{w m_w} h^\dagger_{a m} |0_{\rm core}\rangle ,
  \f]
  with \f$ \eta_{vw} = 1/\sqrt{2} \f$ if v = w (then \f$ J_{12} \f$ is
  even). There is one CSF for each \f$ J_{12} \f$ allowed by the triangle
  rule with \f$ j_a \f$ and J. The parity is the product of the parities of
  the three orbitals. The hole is not antisymmetrised with the particles (the
  hole orbital is in the core, the particles are excited), so no parent
  expansion is needed (cf. @ref CSF3).

  Counts one valence electron (two particles, one hole), and belongs to the
  same (J, parity) blocks as the one-electron CSFs (@ref CSF1): e.g., for Tl
  with \f$ 6s^2 \f$ in the core, \f$ |6p\rangle \f$ and
  \f$ |6p^2\, 6s^{-1}\rangle \f$.

  Particle-hole CI(+MBPT): the hole shell is the top shell of the
  Hartree-Fock core, and the hole's orbital energy and Sigma_1 enter with the
  opposite sign.
   - J. C. Berengut, Phys. Rev. A 94, 012502 (2016)

  @note The particle pair is stored sorted, as in @ref CSF2; this fixes the
        phase convention of all matrix elements.
*/
class CSF2p1h {
public:
  //! The two particle orbitals; always sorted
  std::array<DiracSpinor::Index, 2> states;
  //! The hole orbital (a core orbital)
  DiracSpinor::Index hole;

private:
  int m_twoJ12;
  int m_parity;

public:
  /*!
    @param v,w     The two particle orbitals (nk_index; any order).
    @param a       The hole orbital (nk_index).
    @param twoJ12  Twice the angular momentum of the particle pair, 2*J12.
  */
  CSF2p1h(DiracSpinor::Index v, DiracSpinor::Index w, DiracSpinor::Index a,
          int twoJ12);

  //! Index (nk_index) of the ith particle orbital (i = 0 or 1)
  DiracSpinor::Index state(std::size_t i) const;

  //! Twice the angular momentum of the particle pair, 2*J12
  int twoJ12() const;

  /*!
    @brief Number of orbitals by which the determinants of two CSFs differ:
    the particles of A not in B, plus the holes of B not in A (0 to 3).
    @details
    Matrix elements of one- and two-body operators vanish if more than two
    orbitals differ.
  */
  static int num_different(const CSF2p1h &A, const CSF2p1h &B);

  //! Parity of the CSF, +/-1
  int parity() const { return m_parity; }

  //! Number of electrons in the CSF: particles minus holes
  std::size_t num_electrons() const { return 1; }
  //! Number of particles (excited orbitals) in the CSF
  std::size_t num_particles() const { return 2; }
  //! Number of holes (core orbitals removed) in the CSF
  std::size_t num_holes() const { return 1; }

  //! Configuration as a string, e.g., "6p^2 6s^-1", in relativistic or
  //! non-rel form
  std::string config(bool relativistic = false) const;

  friend bool operator==(const CSF2p1h &A, const CSF2p1h &B);
  friend bool operator!=(const CSF2p1h &A, const CSF2p1h &B);
};

//------------------------------------------------------------------------------
/*!
  @brief Configuration state function (CSF): the type the CI code works with.

  @details
  A CSF is an antisymmetrised N-electron basis state with definite total
  angular momentum (\f$ J^2 \f$, \f$ J_z \f$) and parity \f$\pi\f$.
  Each CSF is built from a combination of \f$ N \f$-electron Slater determinants,
  which are themselves built from N single-particle relativistic orbitals.
  Each Slater determinant has a specific set of single-particle angular
  momentum projections \f$ j_{i,z} = m_i \f$, with \f$ \sum_i m_i = J_z \f$:

  \f[
    |I; J J_z\rangle = \sum_{\{m\}} d_{\{m\}} \,
    |n_1\kappa_1 m_1, n_2\kappa_2 m_2, \ldots, n_N\kappa_N m_N\rangle,
  \f]

  where the kets on the right are (normalised) Slater determinants, and

  \f[
    \sum_{\{m\}} d_{\{m\}}^2 = 1.
  \f]

  The parity is the product from the orbitals,

  \f[
    \pi = \prod_i (-1)^{l_i}.
  \f]

  The determinants are not formed explicitly: the coefficients
  \f$ d_{\{m\}} \f$ are fixed by angular momentum coupling, so matrix elements
  are evaluated directly in the coupled basis.

  Holds one of @ref CSF1 (one electron), @ref CSF2 (two electrons),
  @ref CSF3 (three electrons) or @ref CSF2p1h (two particles and one hole),
  which carry the orbitals and the coupling; the matrix elements are
  evaluated by the functions of each kind (see, e.g., Hab()), and the
  functions taking CSF dispatch to them. A list of CSFs (e.g., of one J/Pi,
  see @ref PsiJPi) may hold CSFs of different kinds with the same number of
  electrons (particles minus holes): one-electron and two-particle one-hole
  CSFs share a block. Other mixed pairs are not implemented (they assert).
*/
class CSF {
  std::variant<CSF1, CSF2, CSF3, CSF2p1h> m_csf;

public:
  CSF(CSF1 csf) : m_csf(std::move(csf)) {}
  CSF(CSF2 csf) : m_csf(std::move(csf)) {}
  CSF(CSF3 csf) : m_csf(std::move(csf)) {}
  CSF(CSF2p1h csf) : m_csf(std::move(csf)) {}

  //! Parity of the CSF, +/-1
  int parity() const;

  //! Number of electrons in the CSF: particles minus holes
  std::size_t num_electrons() const;
  //! Number of particles (excited orbitals) in the CSF
  std::size_t num_particles() const;
  //! Number of holes (core orbitals removed) in the CSF
  std::size_t num_holes() const;

  //! Single-particle configuration as a string, in relativistic or non-rel form
  std::string config(bool relativistic = false) const;

  //! The particle orbitals (nk_index), sorted
  std::vector<DiracSpinor::Index> orbitals() const;

  //! The hole orbitals (nk_index); empty if there are none
  std::vector<DiracSpinor::Index> holes() const;

  //! Index (nk_index) of the ith particle orbital
  DiracSpinor::Index state(std::size_t i) const;

  //! The held CSF1, CSF2, CSF3 or CSF2p1h, for std::visit
  const std::variant<CSF1, CSF2, CSF3, CSF2p1h> &variant() const {
    return m_csf;
  }

  //! Equal if of the same kind, and equal as that kind
  friend bool operator==(const CSF &A, const CSF &B);
  friend bool operator!=(const CSF &A, const CSF &B);
};

/*!
  @brief Number of orbitals by which the determinants of two CSFs differ,
  counting repeated orbitals.
  @details
  The particles of A that are not in B, plus the holes of B that are not in
  A: the orbitals of the full determinant (core, less the holes, plus the
  particles) of A that are absent from that of B. Matrix elements of one-
  and two-body operators vanish if more than two orbitals differ. See
  CSF1::num_different() etc.; for a one-electron and a two-particle one-hole
  CSF, one if the electron is one of the two particles (the other replaces
  the hole), else two. Other mixed kinds assert (see @ref CSF).
*/
int num_different(const CSF &A, const CSF &B);

/*!
  @brief Forms all one-electron CSFs with given total J and parity.
  @details
  These are the orbitals of @p cisp_basis with j = J and the given parity (a
  single kappa), in the order of the basis.

  @param twoJ         Twice the total angular momentum 2J (odd).
  @param parity       Total parity: +1 (even) or -1 (odd).
  @param cisp_basis   Single-particle basis from which CSFs are constructed.
  @param config_space Configurations to include; all by default.
  @return List of all valid one-electron CSFs for the given J and parity.
*/
std::vector<CSF> form_CSF1s(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis,
                            const ConfigurationSpace &config_space = {});

//------------------------------------------------------------------------------
/*!
  @brief Forms all two-electron CSFs with given total J and parity.
  @details
  Iterates over all pairs of single-particle states in @p cisp_basis and
  retains those whose angular momenta can be coupled to total \f$ J = \f$
  @p twoJ /2 and whose combined parity equals @p parity.  Duplicate pairs are
  excluded by construction.

  @param twoJ         Twice the total angular momentum 2J.
  @param parity       Total parity: +1 (even) or -1 (odd).
  @param cisp_basis   Single-particle basis from which CSFs are constructed.
  @param config_space Configurations to include; all by default.
  @return Sorted list of all valid two-electron CSFs for the given J and parity.
*/
std::vector<CSF> form_CSF2s(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis,
                            const ConfigurationSpace &config_space = {});

/*!
  @brief Forms all three-electron CSFs with given total J and parity.
  @details
  Iterates over all distinct triples (with repetition) of single-particle
  states in @p cisp_basis with total parity @p parity, and forms the CSFs of
  each; see @ref form_CSF3s(std::array<DiracSpinor::Index, 3>, int).

  @param twoJ         Twice the total angular momentum 2J (odd).
  @param parity       Total parity: +1 (even) or -1 (odd).
  @param cisp_basis   Single-particle basis from which CSFs are constructed.
  @param config_space Configurations to include; all by default.
  @return All three-electron CSFs for the given J and parity.
*/
std::vector<CSF> form_CSF3s(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis,
                            const ConfigurationSpace &config_space = {});

/*!
  @brief Forms all two-particle one-hole CSFs with given total J and parity.
  @details
  For each hole orbital of @p holes, and each pair of orbitals of
  @p particle_basis (with repetition) of the right parity, one CSF for each
  pair angular momentum J12 allowed by the triangle rule with the hole and
  J (J12 even for two electrons in the same orbital); see @ref CSF2p1h.

  @param twoJ            Twice the total angular momentum 2J (odd).
  @param parity          Total parity: +1 (even) or -1 (odd).
  @param particle_basis  Basis of the two particles.
  @param holes           Core orbitals that may hold the hole.
  @param config_space    Configurations to include; all by default.
  @return All two-particle one-hole CSFs for the given J and parity.
*/
std::vector<CSF> form_CSF2p1hs(int twoJ, int parity,
                               const std::vector<DiracSpinor> &particle_basis,
                               const std::vector<DiracSpinor> &holes,
                               const ConfigurationSpace &config_space = {});

/*!
  @brief Forms all CSFs of @p num_electrons valence electrons with given
  total J and parity.
  @details
  Calls form_CSF1s(), form_CSF2s() or form_CSF3s(); then, if @p holes is not
  empty, appends the hole CSFs of form_CSF2p1hs() (one valence electron
  only: two particles and one hole).

  @param twoJ          Twice the total angular momentum 2J.
  @param parity        Total parity: +1 (even) or -1 (odd).
  @param cisp_basis    Single-particle basis from which CSFs are constructed.
  @param num_electrons Number of valence electrons: 1, 2 or 3.
  @param config_space  Configurations to include; all by default.
  @param holes         Core orbitals that may hold a hole; none by default.
  @param holes_basis   Basis of the particles of the hole CSFs; @p cisp_basis
                       if empty.
  @return All CSFs for the given J and parity.
*/
std::vector<CSF> form_CSFs(int twoJ, int parity,
                           const std::vector<DiracSpinor> &cisp_basis,
                           int num_electrons,
                           const ConfigurationSpace &config_space = {},
                           const std::vector<DiracSpinor> &holes = {},
                           const std::vector<DiracSpinor> &holes_basis = {});

/*!
  @brief Number of pairs of CSFs that can interact: those whose orbitals
  differ by at most two.
  @details
  Counts the pairs (A, B), with B no later than A in @p CSFs (A = B
  included): the elements of the lower triangle of the CI matrix that are not
  zero by the orbitals alone. So an upper bound on the elements stored in the
  sparse CI matrix (construct_Hci()), some of which are zero by angular
  selection. For one or two electrons every pair can interact.

  CSFs of the same configuration (orbitals) are counted together, so the cost
  scales with the square of the number of configurations.

  @param CSFs  List of CSFs.
  @return Number of pairs that can interact.
*/
std::size_t num_interacting_pairs(const std::vector<CSF> &CSFs);

//------------------------------------------------------------------------------

//! The CSFs of one configuration, and the CSFs (of a second list) whose
//! configurations differ from it by at most one orbital
struct OneBodyConnected {
  //! Indices of the CSFs of the configuration
  std::vector<std::size_t> config_CSFs{};
  //! Indices of the CSFs, in the second list, a one-body operator can connect
  //! them to
  std::vector<std::size_t> connected_CSFs{};
};

/*!
  @brief Groups @p CSFAs by configuration, and finds the CSFs of @p CSFBs a
  one-body operator can connect to each.
  @details
  Matrix elements of a one-body operator vanish between CSFs whose
  determinants differ by more than one orbital (see num_different()).
  Configurations within one orbital share all but one of their orbitals:
  @p CSFBs is indexed by its configurations with one orbital (particle or
  hole) removed, so each configuration of @p CSFAs finds its partners
  directly. A configuration with one more particle and one more hole than
  another (a hole excitation of it) is connected to it if it contains it:
  @p CSFBs is also indexed by its configurations as they are, and with one
  particle and one hole removed. The cost scales with the number of
  connected pairs, rather than with @p CSFAs .size() x @p CSFBs .size().

  @param CSFAs  CSFs to be grouped.
  @param CSFBs  CSFs to be connected to them.
  @return One entry per configuration of @p CSFAs.
*/
std::vector<OneBodyConnected> one_body_connected(const std::vector<CSF> &CSFAs,
                                                 const std::vector<CSF> &CSFBs);

//------------------------------------------------------------------------------

/*!
  @brief Matrix element of a two-body operator between coupled pairs of
  orbitals: (a, b, c, d, 2J) -> \f$ \langle (ab)J | O | (cd)J \rangle \f$.
  @details
  Used for the antisymmetrised particle-particle pair elements of the
  Coulomb, Breit and Sigma_2 interactions (e.g., CSF2_Coulomb() bound to a
  Q^k table; see coulomb_pair_element()), and for the elements between
  particle-hole pairs (see @ref particle_hole_sum).
*/
using PairElement =
  std::function<double(DiracSpinor::Index, DiracSpinor::Index,
                       DiracSpinor::Index, DiracSpinor::Index, int)>;

/*!
  @brief Particle-hole part of a two-body operator between two hole CSFs:
  the sum of the interactions of each particle with the hole.
  @details
  The pair of \f$ |(vw)J_{12}, \tilde a; J\rangle \f$ (see @ref CSF2p1h) is
  recoupled to a particle-hole pair,
  \f[
    |(vw)J_{12}, \tilde a; J\rangle = \sum_K R_K(vwa) \,
    |v, (w\tilde a)K; J\rangle ,
  \f]
  \f[
    R_K(vwa) = (-1)^{j_v + j_w + j_a + J} \sqrt{[J_{12}][K]}
    \begin{Bmatrix} j_v & j_w & J_{12} \\ j_a & J & K \end{Bmatrix} ,
  \f]
  and the operator acts on the particle-hole pair with the other particle as
  spectator (which must then coincide between bra and ket). Each of the two
  particles of the bra and of the ket may be the spectator: four terms, the
  swapped ones carrying the pair exchange phase
  \f$ s_{vw} = (-1)^{j_v + j_w - J_{12} + 1} \f$. With the bra
  \f$ |(xy)J_{12}', \tilde b; J\rangle \f$,
  \f[
    \eta_{vw}\eta_{xy} \sum_K \Big[ \delta_{xv} R_K(vwa) R_K(xyb)
    F_K(y\tilde b; w\tilde a) + s_{vw} \delta_{xw} R_K(wva) R_K(xyb)
    F_K(y\tilde b; v\tilde a) + s_{xy} \delta_{yv} R_K(vwa) R_K(yxb)
    F_K(x\tilde b; w\tilde a) + s_{vw} s_{xy} \delta_{yw} R_K(wva)
    R_K(yxb) F_K(x\tilde b; v\tilde a) \Big] ,
  \f]
  where \f$ F_K(y\tilde b; w\tilde a) \f$ = @p ph_element (y, b, w, a, 2K)
  is the element of the operator between the particle-hole pairs
  \f$ |(w\tilde a)K\rangle \f$ (coupled particle first). For a two-body
  interaction (Coulomb, Breit, Sigma_2) it follows from the particle-particle
  pair elements by recoupling (see particle_hole_element()); for
  \f$ L^2 \f$ and \f$ S^2 \f$, under which the hole behaves as a particle,
  it is the plain two-body element (see @ref L2S2_2p1h).

  @param A          The bra CSF (particles x, y; hole b).
  @param B          The ket CSF (particles v, w; hole a).
  @param twoJ       Twice the total angular momentum 2J of both.
  @param ph_element Elements of the operator between particle-hole pairs,
                    (y, b, w, a, 2K) -> \f$ F_K(y\tilde b; w\tilde a) \f$.
  @return The particle-hole sum above.
*/
double particle_hole_sum(const CSF2p1h &A, const CSF2p1h &B, int twoJ,
                         const PairElement &ph_element);

//------------------------------------------------------------------------------

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
  @brief Matrix elements of L^2 and S^2 between two two-electron CSFs of the
  same J (non-relativistic limit).
  @details
  Both are scalars in the two-electron space, and (with the radial orbitals
  of \f$ j = l \pm 1/2 \f$ treated as identical, see @ref LS_amplitude)
  connect only CSFs of the same non-relativistic configuration, in the same
  stored order. Then
  \f[
    \langle x'y'; J | L^2 | xy; J \rangle
    = \sum_{LS} A_{x'y'}(L,S) A_{xy}(L,S) \, L(L+1),
  \f]
  and likewise for \f$ S^2 \f$ with \f$ S(S+1) \f$.

  @param x,y    Orbitals (nk_index, sorted) of the ket CSF.
  @param xp,yp  Orbitals of the bra CSF.
  @param twoJ   Twice the total angular momentum 2J of both CSFs.
  @return Pair {<L^2>, <S^2>}; zero if the configurations differ.
*/
std::pair<double, double> L2S2_pair(DiracSpinor::Index x, DiracSpinor::Index y,
                                    DiracSpinor::Index xp,
                                    DiracSpinor::Index yp, int twoJ);

/*!
  @brief Matrix elements of 2 l_1.l_2 and 2 s_1.s_2 between particle-hole
  pairs, \f$ \langle (y\tilde b)K | 2 l_1 \cdot l_2 | (w\tilde a)K \rangle
  \f$ (non-relativistic limit).
  @details
  Under L and S the hole behaves as a particle (see @ref L2S2_2p1h), so this
  is the two-particle element in product coupling: a particle and a hole are
  never in the same shell, so there is no antisymmetrisation between them.
  Then \f$ 2 l_1 \cdot l_2 = (l_1 + l_2)^2 - l_1^2 - l_2^2 \f$, with
  \f$ (l_1 + l_2)^2 \f$ from @ref L2S2_pair and the one-body terms
  \f$ \delta_{yw}\delta_{ab} [l_w(l_w+1) + l_a(l_a+1)] \f$; likewise for
  the spins, with \f$ s(s+1) = 3/4 \f$.

  @param y,b   Particle and hole of the bra.
  @param w,a   Particle and hole of the ket.
  @param twoK  Twice the angular momentum 2K of the particle-hole pairs.
  @return Pair {<2 l_1.l_2>, <2 s_1.s_2>}; zero if the configurations differ.
*/
std::pair<double, double> L2S2_particle_hole(DiracSpinor::Index y,
                                             DiracSpinor::Index b,
                                             DiracSpinor::Index w,
                                             DiracSpinor::Index a, int twoK);

/*!
  @brief Matrix elements of L^2 and S^2 between two hole CSFs of the same J
  (non-relativistic limit).
  @details
  The closed core has L = S = 0, so (with the radial orbitals of
  \f$ j = l \pm 1/2 \f$ treated as identical, and no mixing between core
  and excited orbitals) the orbital and spin angular momenta act on a hole
  exactly as on a particle: the one-body operators L and S conserve the
  numbers of particles and holes, and their reduced matrix elements between
  hole states equal those between the orbitals. Then, for
  \f$ |(vw)J_{12}, \tilde a; J\rangle \f$ (see @ref CSF2p1h),
  \f$ L^2 = (l_v + l_w + l_a)^2 \f$ is evaluated as for three
  distinguishable particles:
  - the pair, \f$ (l_v + l_w)^2 \f$, with the hole as spectator:
    \f$ \delta_{ab} \delta_{J_{12} J_{12}'} \f$ times @ref L2S2_pair;
  - the hole, \f$ l_a^2 = l_a(l_a+1) \f$, with the pair as spectator;
  - the cross terms, \f$ 2 l_v \cdot l_a + 2 l_w \cdot l_a \f$: the sum
    over the particles of an operator acting on a particle-hole pair,
    evaluated by @ref particle_hole_sum with the elements of
    @ref L2S2_particle_hole.

  Likewise for \f$ S^2 \f$, with \f$ s(s+1) = 3/4 \f$.

  @param A,B   The two CSFs.
  @param twoJ  Twice the total angular momentum 2J of both CSFs.
  @return Pair {<A|L^2|B>, <A|S^2|B>}; zero if the configurations differ.
*/
std::pair<double, double> L2S2_2p1h(const CSF2p1h &A, const CSF2p1h &B,
                                    int twoJ);

/*!
  @brief Matrix elements of L^2 and S^2 between two CSFs of the same J
  (non-relativistic limit).
  @details
  Both are scalars, and (with the radial orbitals of \f$ j = l \pm 1/2 \f$
  treated as identical, see @ref LS_amplitude) connect only CSFs of the same
  non-relativistic configuration, in the same stored order.

  - One electron: each orbital has definite l, and s = 1/2, so both are
    diagonal, \f$ l(l+1) \f$ and \f$ 3/4 \f$.
  - Two electrons: @ref L2S2_pair.
  - Three electrons: written as sums of pair operators plus one-body terms,
  \f[
    L^2 = \sum_{i<j} L_{ij}^2 - \sum_i l_i(l_i+1),
  \f]
  \f[
    S^2 = \sum_{i<j} S_{ij}^2 - \tfrac{9}{4},
  \f]
    the pair terms reduce, by the expansion over parent states (see
    @ref CSF3), to the two-electron matrix elements of @ref L2S2_pair,
    between parents with the same third orbital; the one-body terms are
    diagonal.
  - Hole CSFs (@ref CSF2p1h): @ref L2S2_2p1h. Between a one-electron CSF
    and a hole CSF: zero (L and S conserve the numbers of particles and
    holes).

  Both CSFs must be of the same kind (see @ref CSF), or one a one-electron
  and the other a hole CSF.

  @param A,B   The two CSFs.
  @param twoJ  Twice the total angular momentum 2J of both CSFs.
  @return Pair {<A|L^2|B>, <A|S^2|B>}; zero if the configurations differ.
*/
std::pair<double, double> L2S2_element(const CSF &A, const CSF &B, int twoJ);

/*!
  @brief Expectation values of L^2 and S^2 for a CI state.
  @details
  \f[
    \langle L^2 \rangle = \sum_{IJ} c_I c_J \langle I | L^2 | J \rangle ,
  \f]
  and likewise for \f$ S^2 \f$, with the matrix elements of
  @ref L2S2_element, which connect only CSFs of the same non-relativistic
  configuration.

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
                 const std::vector<CSF> &csfs, int twoJ);

//------------------------------------------------------------------------------

/*!
  @brief Identifies one CI level: its (J, parity), and which solution.
  @details
  The standard text form is `J{+,-}:index`, e.g., `2+:3` is the fourth solution
  (index counts from zero) of the J=2, even parity, CI. 
  The index may be omitted, in which case it is zero:
  `0+` is the lowest even-parity J=0 solution. Half-integer J is written as a
  fraction: `3/2-:1`. See @ref parse_level and @ref to_string.
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
  Accepts `2+:3` (standard) and `e2:3`; the index is optional. Half-integer J
  is given as a fraction over 2, or a decimal: `3/2-:1`, `1.5-:1`, `o3/2:1`.
  Surrounding whitespace is ignored.

  @param str  Text form, e.g., `2+:3`, `e2:3`, `0+`, `1/2-`.
  @return The level; empty if @p str is not a valid level reference.
*/
[[nodiscard]] std::optional<Level> parse_level(std::string_view str);

//! Text form of a CI level reference, e.g., "2+:3", "1/2-:0"; see @ref Level
[[nodiscard]] std::string to_string(const Level &level);

//! J as text, given 2J: "2" for even 2J, "3/2" for odd
[[nodiscard]] std::string J_string(int twoJ);

/*!
  @brief Parses a J value given as text: an integer, a fraction over 2, or a
  decimal.
  @details
  Accepts `2`, `1/2`, `3/2`, `0.5`, `1.5`, `2.0`; a decimal must be an integer
  or half-integer. Surrounding whitespace is ignored.
  @return 2J; empty if @p str is not a valid J.
*/
[[nodiscard]] std::optional<int> parse_twoJ(std::string_view str);

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

  Construction builds the CSF basis (form_CSFs()) but does not solve the
  eigenvalue problem; call solve() separately after constructing the CI
  Hamiltonian matrix.
  Configuration labels are not set automatically -- call update_config_info()
  for each solution after solving.

  The CSFs (see @ref CSF) are those of num_electrons() valence electrons,
  plus the hole CSFs if hole orbitals were given.
*/
class PsiJPi {

  int m_twoj{-1};
  int m_pi{0};
  // Number of valence electrons (particles minus holes) the CSFs were formed
  // with
  int m_num_electrons{0};

  // Number of solutions stored:
  std::size_t m_num_solutions{0};
  std::vector<CSF> m_CSFs{};
  // Energy, and CI expansion coeficients
  std::pair<LinAlg::Vector<double>, LinAlg::Matrix<double>> m_Solution{};
  std::vector<ConfigInfo> m_Info{};

public:
  /*!
    @brief Constructs the CSF basis for the given J and parity; does not solve.
    @details
    Calls form_CSFs() to build the list of CSFs. The eigenvalue problem is
    not solved until solve() is called with the CI Hamiltonian.

    @param twoJ           Twice the total angular momentum 2J.
    @param pi             Total parity: +1 or -1.
    @param cisp_basis     Single-particle basis used to construct the CSFs.
    @param num_electrons  Number of valence electrons: 1, 2 or 3.
    @param config_space   Configurations to include; all by default.
    @param holes          Core orbitals that may hold a hole (see
                          form_CSFs()); none by default.
    @param holes_basis    Basis of the particles of the hole CSFs;
                          @p cisp_basis if empty.
  */
  PsiJPi(int twoJ, int pi, const std::vector<DiracSpinor> &cisp_basis,
         int num_electrons = 2, const ConfigurationSpace &config_space = {},
         const std::vector<DiracSpinor> &holes = {},
         const std::vector<DiracSpinor> &holes_basis = {});

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
    @brief Solves the CI eigenvalue problem iteratively, for a sparse
    Hamiltonian matrix.
    @details
    Finds the lowest @p num_solutions eigenpairs with LinAlg::Davidson, and
    stores them, as solve() does. For large J/Pi blocks, where the dense
    matrix does not fit in memory (or takes too long to diagonalise) and only
    a few solutions are needed.

    If @p all_below is set, all eigenpairs below that energy are found
    instead: their number is not known in advance, so the lowest few are
    found, and the number is doubled (solving again) while the highest found
    is still below the threshold.

    Davidson finds fewer solutions than CSFs: if num_CSFs() or more are
    requested (or all), num_CSFs() - 1 are found. The Davidson vectors grow
    with the number of solutions: beyond ~num_CSFs()/20 it is slower, and
    takes more memory, than the dense solve(), and a warning is printed. A
    single CSF is its own solution.

    @param Hci           CI Hamiltonian matrix in the CSF basis (sparse; see
                         construct_Hci()).
    @param num_solutions Number of lowest solutions to find (all if <= 0).
                         Ignored if @p all_below is set.
    @param all_below     If set, find all solutions below this energy (cm^-1).
    @param outstream     For the warning when many solutions are requested.
    @return Convergence of the iterative solution (of the last Davidson
            solve, if repeated); {0 iterations, 0 residual, converged} for a
            single CSF.
  */
  LinAlg::DavidsonConvergence
  solve(const LinAlg::SparseSymmetricMatrix<double> &Hci, int num_solutions,
        std::optional<double> all_below = {},
        std::ostream &outstream = std::cout);

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

  /*!
    @brief Replaces the ith solution (energy and expansion coefficients).
    @details
    For a level solved separately from the others: e.g., with its own reference
    energy E0 (see @ref run_CI). The solutions are then no longer eigenstates
    of a single Hamiltonian, so are not orthogonal.

    @param i      Index of an existing solution.
    @param energy Its new energy (atomic units).
    @param coefs  Its new CI expansion coefficients; one per CSF.
  */
  void replace_solution(std::size_t i, double energy,
                        const LinAlg::View<const double> &coefs);

  //! Set configuration info for the ith solution (must be called manually after solve())
  void update_config_info(std::size_t i, const ConfigInfo &info);

  //! Number of valence electrons (particles minus holes) the CSFs were
  //! formed with (1, 2 or 3)
  int num_electrons() const;

  //! True if any of the CSFs has a hole
  bool has_holes() const;

  //! Number of CSFs spanning this J/Pi
  std::size_t num_CSFs() const;

  //! The ith CSF
  const CI::CSF &csf(std::size_t i) const;

  //! Full list of CSFs spanning this J/Pi
  const std::vector<CSF> &CSFs() const;

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
    A single file holds multiple blocks, one per (twoJ, parity) pair, after a
    short format tag. Each block is self-describing:
    (twoJ, pi, num_csfs, num_solutions, E[num_solutions],
    M[num_solutions x num_csfs]), with row i of M the CSF coefficients of
    solution i. Only the solutions found are stored, so the file size does
    not grow as num_csfs^2 when only a few solutions are found.

    The CSF basis (m_CSFs) is not touched; it must be constructed via the normal
    constructor before calling this function. The basis (and hence num_csfs) must
    be consistent with the file on write; a mismatch causes failure.

    ConfigInfo is not stored -- call update_config_info() after reading if needed.

    On read: checks the format tag (files without it are not read), scans
    blocks until the matching (twoJ, pi) is found, verifies num_csfs, reads
    num_solutions, E, and M; returns false if the block is not found or
    num_csfs mismatches.

    On write: a new (twoJ, pi) block is appended. If the block exists already
    (its size may differ), the file is re-written: the other blocks, then
    this one. A file without the format tag is replaced.

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
