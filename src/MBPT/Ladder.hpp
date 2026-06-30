#pragma once
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"

namespace MBPT {

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
  @return \f$ L^k_{mnij} \f$
  @note To evaluate at a fixed external energy (for a correlation potential),
        pass the external orbital with its energy set accordingly: the energy
        only enters the denominators, never the integral lookups.
*/
double Lkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited, bool include_L4,
              const Angular::SixJTable &SJ,
              const Coulomb::LkTable *const Lk = nullptr);

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
  @return \f$ L1^k_{mnij} \f$
*/
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr);

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
  @return \f$ L2^k_{mnij} \f$
*/
double L2(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr);

/*!
  @brief Exchange partner of L2; equals L2 with m,n and i,j swapped.
  @details
  \f[ L3^k_{mnij} = L2^k_{nmji} \f]
*/
inline double L3(int k, const DiracSpinor &m, const DiracSpinor &n,
                 const DiracSpinor &i, const DiracSpinor &j,
                 const Coulomb::QkTable &qk,
                 const std::vector<DiracSpinor> &core,
                 const std::vector<DiracSpinor> &excited,
                 const Angular::SixJTable &SJ,
                 const Coulomb::LkTable *const Lk = nullptr) {
  return L2(k, n, m, j, i, qk, core, excited, SJ, Lk);
}

/*!
  @brief Core--core (hole--hole) ladder diagram L4.
  @details
  Intermediate states run over core orbitals only, making this the
  hole--hole counterpart of the particle--particle diagram L1.
  Typically small and omitted by default; enable via @p include_L4 in
  Lkmnij().

  @param k   Multipole rank
  @param m,n Excited (particle) orbitals
  @param i,j Core (hole) or valence orbitals
  @param qk  Coulomb \f$ Q^k \f$ integral table
  @param core Core orbitals
  @param SJ  6j symbol table
  @param Lk  Ladder table from previous iteration (nullptr on first)
  @return \f$ L4^k_{mnij} \f$
*/
double L4(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr);

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
  @brief Ladder-diagram correction to the correlation potential, Sigma_L(en_v).
  @details
  Forms the ladder correlation potential by projecting the discrete ladder
  integrals onto the excited basis states of @p kappa_v. The exchange is folded
  into the Coulomb vertex via \f$ W = Q + P \f$ (as in de_valence_w):
  \f[
    \Sigma_L = \sum_{i,amn,k} |W^k_{\cdot amn}\rangle\,
      \frac{L^k_{mn,i,a}}{[k][j_v]\,(\epsilon_v+\epsilon_a-\epsilon_m-\epsilon_n)}\,
      \langle i|
      + \sum_{i,nab,k} |W^k_{\cdot nab}\rangle\,
      \frac{L^k_{i,n,a,b}}{[k][j_v]\,(\epsilon_v+\epsilon_n-\epsilon_a-\epsilon_b)}\,
      \langle i| ,
  \f]
  (particle-particle (a+b) and particle-hole (c+d) diagrams). The bra index
  \f$ i \f$ runs over the excited basis states of @p kappa_v (approximating
  completeness). The ladder integrals are computed on-the-fly via Lkmnij()
  evaluated at the fixed external energy @p en_v (applied by setting the energy
  of the external orbital, since the energy enters only the denominators).

  The sub-grid (@p r0, @p rmax, @p stride) defaults match Wavefunction::formSigma.

  @param kappa_v   Valence kappa
  @param en_v      Energy at which to evaluate Sigma_L
  @param core      Core (hole) orbitals
  @param excited   Excited orbitals (also supplies the projection basis)
  @param qk        Converged Coulomb \f$ Q^k \f$ table
  @param lk        Internal-rung ladder table, forwarded to Lkmnij: nullptr for
                   L(Q,Q)=L^(1); a converged table (its fixed point) for the
                   full ladder
  @param sjt       6j symbol table
  @param include_L4 Include core--core diagram in on-the-fly Lkmnij
  @param r0,rmax,stride  Sub-grid parameters
  @param include_G Include the lower (g) component of Sigma_L
  @return Sigma_L as a coordinate-space GMatrix
*/
GMatrix Sigma_ladder(int kappa_v, double en_v,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited,
                     const std::vector<DiracSpinor> &basis,
                     const Coulomb::QkTable &qk, const Coulomb::LkTable *lk,
                     const Angular::SixJTable &sjt, bool include_L4 = false,
                     double r0 = 1.0e-4, double rmax = 30.0,
                     std::size_t stride = 4, bool include_G = false);

// template implementations:
#include "Ladder.ipp"

} // namespace MBPT
