#pragma once
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <optional>
#include <string>
#include <utility>

namespace MBPT {

//! Returns min and max k (multipolarity) allowed for ladder integral
//! L^k_abcd. Triangle rules {a,c,k} and {b,d,k} apply, plus the combined
//! parity rule (l_a+l_b+l_c+l_d even). Unlike Q^k, there is no
//! individual-pair parity rule linking k to the orbital parities, so k takes
//! both parities: NOT safe to call k+=2.
inline std::pair<int, int> k_minmax_L(const DiracSpinor &a,
                                      const DiracSpinor &b,
                                      const DiracSpinor &c,
                                      const DiracSpinor &d) {
  // Combined parity rule: each (direct) Coulomb rung shifts the parity of
  // both pair lines together, so only the total parity is constrained
  if ((a.l() + b.l() + c.l() + d.l()) % 2 != 0) {
    return {1, 0};
  }
  const auto [l1, u1] = Coulomb::k_minmax_tj(a.twoj(), c.twoj());
  const auto [l2, u2] = Coulomb::k_minmax_tj(b.twoj(), d.twoj());
  return {std::max(l1, l2), std::min(u1, u2)};
}

//! Returns min and max k (multipolarity) allowed for ladder integral L^k_abcd
//! (kappa version) - see above. NOT safe to call k+=2.
inline std::pair<int, int> k_minmax_L(int kap_a, int kap_b, int kap_c,
                                      int kap_d) {
  if ((Angular::l_k(kap_a) + Angular::l_k(kap_b) + Angular::l_k(kap_c) +
       Angular::l_k(kap_d)) %
        2 !=
      0) {
    return {1, 0};
  }
  const auto [l1, u1] =
    Coulomb::k_minmax_tj(Angular::twoj_k(kap_a), Angular::twoj_k(kap_c));
  const auto [l2, u2] =
    Coulomb::k_minmax_tj(Angular::twoj_k(kap_b), Angular::twoj_k(kap_d));
  return {std::max(l1, l2), std::min(u1, u2)};
}

/*!
  @brief Full ladder integral summed over all diagrams.
  @details
  Computes
  \f[
    L^k_{mnij} = L1^k_{mnij} + L2^k_{mnij} + L3^k_{mnij} [+ L4^k_{mnij}]
  \f]
  where \f$ L3^k_{mnij} = L2^k_{nmji} \f$ and \f$ L4 \f$ involves
  core--core intermediate states. @p Lk points to the ladder table from
  the previous iteration; pass nullptr on the first iteration.

  @param k          Multipole rank
  @param m,n        Excited (particle) orbitals
  @param i,j        Core (hole) or valence orbitals
  @param qk         Coulomb \f$ Q^k \f$ integral table
  @param core       Core orbitals
  @param excited    Excited orbitals
  @param include_L4 Include the core--core diagram L4
  @param SJ         6j symbol table
  @param Lk         Ladder table from previous iteration (nullptr on first)
  @param e_i        Optional: used in place of i.en() in energy denominators
  @param e_m        Optional: used in place of m.en() in energy denominators
  @return \f$ L^k_{mnij} \f$
  @note To evaluate at a fixed external energy (for a correlation potential),
        use @p e_i or @p e_m (for external line in the i or m slot): the energy
        only enters the denominators, never the integral lookups.
*/
double Lkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited, bool include_L4,
              const Angular::SixJTable &SJ,
              const Coulomb::LkTable *const Lk = nullptr,
              std::optional<double> e_i = {}, std::optional<double> e_m = {});

/*!
  @brief Particle--particle ladder diagram L1.
  @details
  \f[
    L1^k_{mnij} = \sum_{rs,ul} A^{kul}_{mnrsij}
      \frac{Q^u_{mnrs}\,(Q+L)^l_{rsij}}{\epsilon_{ij} - \epsilon_{rs}}
  \f]
  with the angular coefficient
  \f[
    A^{kul}_{mnrsij} = (-1)^{m+n+r+s+i+j+1}\,[k]
      \sixj{m}{i}{k}{l}{u}{r}\sixj{n}{j}{k}{l}{u}{s}
  \f]
  Intermediate states \f$ r,s \f$ run over excited orbitals.

  @param k       Multipole rank
  @param m,n     Excited (particle) orbitals
  @param i,j     Core (hole) or valence orbitals
  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param excited Excited orbitals
  @param SJ      6j symbol table
  @param Lk      Ladder table from previous iteration (nullptr on first)
  @param e_i     Optional: used in place of i.en() in energy denominator
  @return \f$ L1^k_{mnij} \f$
*/
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr,
          std::optional<double> e_i = {});

/*!
  @brief Particle--hole ladder diagram L2.
  @details
  \f[
    L2^k_{mnij} = \sum_{rc,ul} (-1)^{k+u+l+1} A^{klu}_{mjcrin}
      \frac{Q^u_{cnir}\,(Q+L)^l_{mrcj}}{\epsilon_{cj} - \epsilon_{mr}}
  \f]
  Intermediate states: \f$ r \f$ runs over excited, \f$ c \f$ over core.
  The diagram \f$ L3 \f$ is the exchange partner \f$ L3^k_{mnij} = L2^k_{nmji} \f$.

  @param k       Multipole rank
  @param m,n     Excited (particle) orbitals
  @param i,j     Core (hole) or valence orbitals
  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param core    Core orbitals
  @param excited Excited orbitals
  @param SJ      6j symbol table
  @param Lk      Ladder table from previous iteration (nullptr on first)
  @param e_j     Optional: used in place of j.en() in energy denominator
  @param e_m     Optional: used in place of m.en() in energy denominator
  @return \f$ L2^k_{mnij} \f$
*/
double L2(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr,
          std::optional<double> e_j = {}, std::optional<double> e_m = {});

/*!
  @brief Exchange partner of L2; equals L2 with m,n and i,j swapped.
  @details
  \f[ L3^k_{mnij} = L2^k_{nmji} \f]
  (@p e_i optional: used in place of i.en() in energy denominator)
*/
inline double
L3(int k, const DiracSpinor &m, const DiracSpinor &n, const DiracSpinor &i,
   const DiracSpinor &j, const Coulomb::QkTable &qk,
   const std::vector<DiracSpinor> &core,
   const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
   const Coulomb::LkTable *const Lk = nullptr, std::optional<double> e_i = {}) {
  return L2(k, n, m, j, i, qk, core, excited, SJ, Lk, e_i);
}

/*!
  @brief Core--core (hole--hole) ladder diagram L4.
  @details
  Intermediate states run over core orbitals only, making this the
  hole--hole counterpart of the particle--particle diagram L1.
  Enable via @p include_L4 in Lkmnij().

  @param k   Multipole rank
  @param m,n Excited (particle) orbitals
  @param i,j Core (hole) or valence orbitals
  @param qk  Coulomb \f$ Q^k \f$ integral table
  @param core Core orbitals
  @param SJ  6j symbol table
  @param Lk  Ladder table from previous iteration (nullptr on first)
  @param e_m Optional: used in place of m.en() in energy denominator
  @return \f$ L4^k_{mnij} \f$
*/
double L4(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr,
          std::optional<double> e_m = {});

/*!
  @brief Fills the ladder integral table for all _new_ index combinations.
  @details
  Iterates over all combinations of excited pairs \f$ (m,n) \f$ and
  orbitals in @p i_orbs, computing \f$ L^k_{mnib} \f$ and storing results
  in @p lk.
  Only calculates new integrals. Only lowest-order.

  @param lk          Output ladder table (written in place)
  @param qk          Coulomb \f$ Q^k \f$ integral table
  @param excited      Excited orbitals
  @param core        Core orbitals
  @param i_orbs      Orbitals for the \f$ i \f$ index
  @param include_L4  Include core-core diagram L4
  @param sjt         6j symbol table
  @param max_k       Maximum multipolarity; -1 uses qk.max_k()
  @param print       Print Qk info to screen
*/
void fill_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &excited,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &i_orbs, bool include_L4,
                  const Angular::SixJTable &sjt, int max_k = -1,
                  bool print = true);

/*!
  @brief Updates the ladder integral table with L(Q,Q) -> L(Q,Q+L)
  @details
  Iterates over all combinations of excited pairs \f$ (m,n) \f$ and
  orbitals in @p i_orbs, computing \f$ L^k_{mnib} \f$ and storing results
  in @p lk. Designed for iterative refinement: pass the previous iteration's
  table as @p lk_prev.

  @note Does not calculate any new integrals - assumes all already present.
  Just updates them (based on iterative rule: L(Q,Q) -> L(Q,Q+L))

  @param lk          Output ladder table (written in place)
  @param qk          Coulomb \f$ Q^k \f$ integral table
  @param excited      Excited orbitals
  @param core        Core orbitals
  @param update_i    Restrict re-iteration to entries whose i index is in this
                     set (b is always core). Empty => update all. Used to
                     converge core (update_i=core) before valence
                     (update_i=valence).
  @param include_L4  Include core--core diagram L4
  @param sjt         6j symbol table
  @param lk_prev     Ladder table from previous iteration
  @param a_damp      Damping factor [0,1) : 0 means no damping
  @param print       Print Qk info to screen
*/
void update_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                    const std::vector<DiracSpinor> &excited,
                    const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &update_i, bool include_L4,
                    const Angular::SixJTable &sjt,
                    const Coulomb::LkTable *const lk_prev, double a_damp,
                    bool print);

/*!
  @brief Second-order (or ladder) correction to the valence energy.
  @details
  Computes the correlation energy shift
  \f[
    \delta\epsilon_v = \sum_{mnc} \frac{Q^k_{vmcn}\,L^k_{mncv}}{\epsilon_v
      + \epsilon_c - \epsilon_m - \epsilon_n}
  \f]
  (schematic). When @p lk holds plain Coulomb integrals the result is the
  MBPT(2) correction; when @p lk holds ladder integrals it is the full
  ladder correction.

  @param v       Valence orbital
  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param lk      \f$ Q^k \f$ or ladder \f$ L^k \f$ integral table
  @param core    Core orbitals
  @param excited Excited orbitals
  @return \f$ \delta\epsilon_v \f$
*/
template <typename Qintegrals, typename QorLintegrals>
double de_valence(const DiracSpinor &v, const Qintegrals &qk,
                  const QorLintegrals &lk, const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited);

/*!
  @brief Ladder (or MBPT2) valence energy, antisymmetrising the FIRST integral.
  @details
  Computes
  \f[
    \delta\epsilon_v = \sum_{amn,k}\frac{W^k_{vamn}\,L^k_{mnva}}{[k]\,[j_v]\,\Delta\epsilon}
      + \text{(c+d)},
  \f]
  with \f$ W = Q + P \f$ the antisymmetrised Coulomb integral (@p qk.W) and the
  ladder @p lk entering only through the direct integral @p lk.Q. Equivalent to
  de_valence() (which instead antisymmetrises the ladder via @p lk.P), since the
  exchange symmetry holds under the full sum. Unlike de_valence(), this works
  only with the Coulomb integrals in the first slot (the ladder lacks the
  required symmetry). No screening/eta is applied.

  @param v       Valence orbital
  @param qk      Coulomb integrals supplying \f$ W = Q+P \f$ (first integral)
  @param lk      Ladder \f$ L^k \f$ integrals (second, direct integral)
  @param core    Core orbitals
  @param excited Excited orbitals
  @param sj      Optional 6j table (speeds up the W exchange sum)
  @return \f$ \delta\epsilon_v \f$
*/
template <typename Qintegrals, typename Lintegrals>
double de_valence_w(const DiracSpinor &v, const Qintegrals &qk,
                    const Lintegrals &lk, const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &excited,
                    const Angular::SixJTable *sj = nullptr);

/*!
  @brief Second-order (or ladder) correction to the core energy.
  @details
  Sums the correlation energy shift over all core orbitals. When @p lk holds
  plain Coulomb integrals the result is the MBPT(2) core correction; when
  @p lk holds ladder integrals it is the ladder correction.

  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param lk      \f$ Q^k \f$ or ladder \f$ L^k \f$ integral table
  @param core    Core orbitals
  @param excited Excited orbitals
  @return Total core correlation energy shift
*/
template <typename Qintegrals, typename QorLintegrals>
double de_core(const Qintegrals &qk, const QorLintegrals &lk,
               const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited);

/*!
  @brief Method used to construct the ladder correlation potential, Sigma_L.
  @details
  - single : project onto the single HF |v> eigenstate [Sigma_ladder()]
  - ladder : project onto states in the ladder basis [Sigma_ladder()]
  - full   : project onto the entire basis; requires extending Qk (slow)
             [Sigma_ladder()]
  - ratio  : no projection; rescale each Sigma(2) term by L/Q
             [Sigma_ladder(), empty projection basis]
  - direct : no projection; open the external line exactly (ladder vertex)
             [Sigma_ladder_direct()]
*/
enum class SigmaLMethod { single, ladder, full, ratio, direct };

//! Converts string (name) to SigmaLMethod enum (case-insensitive); warns and
//! defaults to ladder if unknown
SigmaLMethod parseSigmaLMethod(const std::string &method);
//! Converts SigmaLMethod enum to string (name)
std::string parseSigmaLMethod(SigmaLMethod method);

/*!
  @brief Ladder-diagram correction to the correlation potential, Sigma_L(e_v),
  by projection (or, if @p projection is empty, by the L/Q ratio method).
  @details
  With a non-empty @p projection basis, forms the ladder correlation potential
  by projecting the discrete ladder integrals onto the projection states of
  kappa_v. The exchange is folded into the Coulomb vertex via
  \f$ W = Q + P \f$ (as in de_valence_w):
  \f[
    \Sigma_L = \sum_{i,amn,k} |W^k_{\cdot amn}\rangle\,
      \frac{L^k_{mn,i,a}}{[k][j_v]\,(\epsilon_v+\epsilon_a-\epsilon_m-\epsilon_n)}\,
      \langle i|
      + \sum_{i,nab,k} |W^k_{\cdot nab}\rangle\,
      \frac{L^k_{i,n,a,b}}{[k][j_v]\,(\epsilon_v+\epsilon_n-\epsilon_a-\epsilon_b)}\,
      \langle i| ,
  \f]
  (particle-particle (a+b) and particle-hole (c+d) diagrams). The bra index
  \f$ i \f$ runs over the @p projection states of kappa_v (approximating
  completeness). The ladder integrals are computed on-the-fly via Lkmnij()
  evaluated at the fixed external energy \f$ \epsilon_v \f$ (via the e_i/e_m
  energy overrides, since the energy enters only the denominators). Exception:
  for the valence state itself, the stored table entries are used directly -
  they are already at the correct energy, making single-state projection
  essentially free.

  If @p projection is empty, instead uses the ratio method (following
  V. A. Dzuba, Phys. Rev. A 78, 042502 (2008)): no projection; each term of
  the regular second-order correlation potential (cf. Goldstone::Sigma_both)
  is rescaled by the scalar ratio of ladder to Coulomb integrals:
  \f[
    \Sigma_L = \sum_{amn,k} |Q^k_{\cdot amn}\rangle\,
      \frac{L^k_{mnva}/Q^k_{mnva}}{[k][j_v]\,(\epsilon_v+\epsilon_a-\epsilon_m-\epsilon_n)}\,
      \langle W^k_{\cdot amn}|
      + \sum_{nab,k} |Q^k_{\cdot nab}\rangle\,
      \frac{L^k_{vnab}/Q^k_{vnab}}{[k][j_v]\,(\epsilon_v+\epsilon_n-\epsilon_a-\epsilon_b)}\,
      \langle W^k_{\cdot nab}| .
  \f]
  By construction the diagonal reproduces the ladder energy exactly:
  <v|Sigma_L|v> = de_valence_w(v). The off-diagonal (radial) structure is
  approximate: each term keeps the shape of the corresponding second-order
  term, rescaled by a scalar. All integrals come straight from the stored
  tables (the L^k entries with i = v are already at the valence energy), so
  nothing is computed on-the-fly: much faster than projection.

  The sub-grid (@p r0, @p rmax, @p stride) defaults match Wavefunction::formSigma.

  @param v         Valence state (basis version; must be in the qk/lk tables)
  @param core      Core (hole) orbitals
  @param excited   Excited orbitals
  @param projection Projection basis {|i>}; states of kappa_v are used. If
                   empty, the ratio method is used instead (no projection)
  @param qk        Converged Coulomb \f$ Q^k \f$ table
  @param lk        Converged ladder \f$ L^k \f$ table. For projection, it is
                   the internal-rung table forwarded to Lkmnij (nullptr for
                   L(Q,Q)=L^(1)); for ratio, it supplies the L^k integrals
                   (nullptr gives Sigma_L = 0)
  @param sjt       6j symbol table
  @param include_L4 Include core--core diagram in on-the-fly Lkmnij
                   (projection only; the ratio method never re-computes L)
  @param r0,rmax,stride  Sub-grid parameters
  @param include_G Include the lower (g) component of Sigma_L
  @return Sigma_L as a coordinate-space GMatrix

  @note Ratio method: terms with \f$ Q^k_{mnva} = 0 \f$ (but \f$ L^k \ne 0 \f$)
        cannot be rescaled and are dropped; inherent to the method. Note that
        L^k exists at both parities of k (see k_minmax_L) while Q^k
        exists only at Coulomb parity, so the exchange-only (wrong-parity)
        k channels are always dropped here; the projection and direct methods
        include them.
*/
GMatrix Sigma_ladder(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited,
                     const std::vector<DiracSpinor> &projection,
                     const Coulomb::QkTable &qk, const Coulomb::LkTable *lk,
                     const Angular::SixJTable &sjt, bool include_L4 = false,
                     double r0 = 1.0e-4, double rmax = 30.0,
                     std::size_t stride = 4, bool include_G = false);

/*!
  @brief Vertex (ket) form of the ladder integral L^k_mnia over the external
  index i.
  @details
  Returns the radial spinor \f$ |L^k_{mn \cdot a}\rangle \f$ (kappa of @p v)
  satisfying
  \f[ \langle x|L^k_{mn\cdot a}\rangle = L^k_{mnxa}(\epsilon_v) \f]
  for any \f$ x \f$ with \f$ \kappa_x = \kappa_v \f$, with the external energy
  fixed at \f$ \epsilon_v \f$ (as the e_i override in Lkmnij).

  In each diagram the external line attaches to a single bare Coulomb line;
  that line is opened exactly as a radial function (Qkv_bcd). The exception is
  the L-part of the internal (Q+L) rung in L1 and L3, where the external line
  attaches to a dressed rung: for that piece we set i = v (scalar coefficient
  times |v>), which is exact at x = v and is the same level of treatment the
  scalar table gives it (entries exist only for stored orbitals).

  @note All internal lines are Hartree-Fock basis states (from the tables);
        only the external line is left open. This is the correct structure for
        acting on Brueckner orbitals.
*/
DiracSpinor Lkv_mnia(int k, const DiracSpinor &v, const DiracSpinor &m,
                     const DiracSpinor &n, const DiracSpinor &a,
                     const Coulomb::QkTable &qk, const Coulomb::YkTable &yk,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited, bool include_L4,
                     const Angular::SixJTable &SJ,
                     const Coulomb::LkTable *const Lk = nullptr);

/*!
  @brief Vertex (ket) form of the ladder integral L^k_inab over the external
  index i (the m-slot).
  @details
  Returns the radial spinor \f$ |L^k_{\cdot nab}\rangle \f$ (kappa of @p v)
  satisfying
  \f[ \langle x|L^k_{\cdot nab}\rangle = L^k_{xnab}(\epsilon_v) \f]
  for any \f$ x \f$ with \f$ \kappa_x = \kappa_v \f$, with the external energy
  fixed at \f$ \epsilon_v \f$ (as the e_m override in Lkmnij).

  Mirror of Lkv_mnia() for the particle-hole (c+d) diagrams: here the external
  line sits in the m-slot, so it is L2 and L4 whose internal (Q+L) rung
  contains it (i = v used for the L-part), while L1 and L3 open exactly.
*/
DiracSpinor Lkv_inab(int k, const DiracSpinor &v, const DiracSpinor &n,
                     const DiracSpinor &a, const DiracSpinor &b,
                     const Coulomb::QkTable &qk, const Coulomb::YkTable &yk,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited, bool include_L4,
                     const Angular::SixJTable &SJ,
                     const Coulomb::LkTable *const Lk = nullptr);

/*!
  @brief Ladder correlation potential Sigma_L via the direct (open external
  line) method.
  @details
  Forms the ladder correlation potential with the external line opened
  exactly, rather than projected onto a basis:
  \f[
    \Sigma_L = \sum_{amn,k} |W^k_{\cdot amn}\rangle\,
      \frac{1}{[k][j_v]\,(\epsilon_v+\epsilon_a-\epsilon_m-\epsilon_n)}\,
      \langle L^k_{mn \cdot a}|
      + \sum_{nab,k} |W^k_{\cdot nab}\rangle\,
      \frac{1}{[k][j_v]\,(\epsilon_v+\epsilon_n-\epsilon_a-\epsilon_b)}\,
      \langle L^k_{\cdot nab}| ,
  \f]
  with the bra-side ladder vertices from Lkv_mnia() / Lkv_inab(). All internal
  lines are Hartree-Fock basis states; the external line is exact wherever it
  attaches to a bare Coulomb line (everything at lowest order in L), and is
  taken as |v><v| only for the dressed-rung (internal-L) attachment, which
  enters the energy at 4th order. Consequently Sigma_L acts correctly on
  Brueckner orbitals: (H + Sigma + Sigma_L)|psi_B> = e|psi_B>, rather than
  merely shifting by <v|Sigma_L|v>.

  <v|Sigma_L|v> reproduces the ladder energy de_valence_w (up to iteration
  convergence of the L table). No projection basis, no Qk extension.

  @param v        Valence state (basis version; supplies kappa_v, en_v, and
                  the i=v dressed-rung piece)
  @param core     Core (hole) orbitals
  @param excited  Excited orbitals
  @param qk       Converged Coulomb \f$ Q^k \f$ table
  @param yk       Yk table spanning core+excited (radial vertex functions)
  @param lk       Converged ladder \f$ L^k \f$ table (nullptr: lowest-order L)
  @param sjt      6j symbol table
  @param include_L4 Include core--core diagram L4
  @param r0,rmax,stride  Sub-grid parameters
  @param include_G Include the lower (g) component of Sigma_L
  @return Sigma_L as a coordinate-space GMatrix
*/
GMatrix Sigma_ladder_direct(
  const DiracSpinor &v, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &excited, const Coulomb::QkTable &qk,
  const Coulomb::YkTable &yk, const Coulomb::LkTable *lk,
  const Angular::SixJTable &sjt, bool include_L4 = false, double r0 = 1.0e-4,
  double rmax = 30.0, std::size_t stride = 4, bool include_G = false);

//==============================================================================

//! Ladder correlation potential matrix, Sigma_L, for a single (kappa, n)
//! state, evaluated at energy en.
struct SigmaLData {
  int kappa;
  int n;
  double en;
  GMatrix SL;
};

/*!
  @brief Writes Sigma_L (ladder) matrices to binary file.
  @details
  File contains the full-grid parameters (for checking on read), followed by
  each Sigma_L matrix with its own sub-grid parameters and include_G flag.
  Returns false (and writes nothing) if fname is empty or 'false'.
*/
bool write_SigmaL(const std::string &fname, const std::vector<SigmaLData> &SLs,
                  const Grid &grid);

/*!
  @brief Reads Sigma_L (ladder) matrices from binary file.
  @details
  Returns empty vector if file doesn't exist or on grid mismatch: @p grid must
  match the full radial grid the matrices were calculated on. Each matrix
  carries its own sub-grid parameters and include_G flag (they need not match
  those of the base Sigma).
*/
std::vector<SigmaLData> read_SigmaL(const std::string &fname,
                                    const std::shared_ptr<const Grid> &grid);

// template implementations:
#include "Ladder.ipp"

} // namespace MBPT
