#pragma once
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"

namespace MBPT {

// // Weighted average:
// fk = 0.714, 0.617, 0.832, 0.885, 0.941, 0.970, 0.987, 0.994
// fk = 0.654, 0.514, 0.792, 0.857, 0.930, 0.967, 0.986, 0.994 // w/ hp
// eta= 1.258, 1.511, 1.281, 1.175, 1.114, 1.073, 1.062, 1.066

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
  @param fk         Optional screening factors \f$ f_k \f$
  @return \f$ L^k_{mnij} \f$
*/
double Lkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited, bool include_L4,
              const Angular::SixJTable &SJ,
              const Coulomb::LkTable *const Lk = nullptr,
              const std::vector<double> &fk = {});

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
  @param fk      Optional screening factors \f$ f_k \f$
  @return \f$ L1^k_{mnij} \f$
*/
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr,
          const std::vector<double> &fk = {});

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
  @param fk      Optional screening factors \f$ f_k \f$
  @return \f$ L2^k_{mnij} \f$
*/
double L2(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr,
          const std::vector<double> &fk = {});

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
                 const Coulomb::LkTable *const Lk = nullptr,
                 const std::vector<double> &fk = {}) {
  return L2(k, n, m, j, i, qk, core, excited, SJ, Lk, fk);
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
  @param fk  Optional screening factors \f$ f_k \f$
  @return \f$ L4^k_{mnij} \f$
*/
double L4(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr,
          const std::vector<double> &fk = {});

/*!
  @brief Fills the ladder integral table for all required index combinations.
  @details
  Iterates over all combinations of excited pairs \f$ (m,n) \f$ and
  orbitals in @p i_orbs, computing \f$ L^k_{mnib} \f$ and storing results
  in @p lk. Designed for iterative refinement: pass the previous iteration's
  table as @p lk_prev.

  @param lk          Output ladder table (written in place)
  @param qk          Coulomb \f$ Q^k \f$ integral table
  @param excited      Excited orbitals
  @param core        Core orbitals
  @param i_orbs      Orbitals for the \f$ i \f$ index
  @param include_L4  Include core--core diagram L4
  @param sjt         6j symbol table
  @param lk_prev     Ladder table from previous iteration (nullptr on first)
  @param print_progbar Print a progress bar to stdout
  @param fk          Optional screening factors \f$ f_k \f$
*/
void fill_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &excited,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &i_orbs, bool include_L4,
                  const Angular::SixJTable &sjt,
                  const Coulomb::LkTable *const lk_prev = nullptr,
                  bool print_progbar = true,
                  const std::vector<double> &fk = {});

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
  @param fk      Optional screening factors \f$ f_k \f$
  @param etak    Optional enhancement factors \f$ \eta_k \f$
  @return \f$ \delta\epsilon_v \f$
*/
template <typename Qintegrals, typename QorLintegrals>
double de_valence(const DiracSpinor &v, const Qintegrals &qk,
                  const QorLintegrals &lk, const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const std::vector<double> &fk = {},
                  const std::vector<double> &etak = {});

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

// template implementations:
#include "Ladder.ipp"

} // namespace MBPT
