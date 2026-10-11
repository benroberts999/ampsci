#pragma once
#include "Coulomb/YkTable.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <memory>
#include <vector>

namespace MBPT {

//------------------------------------------------------------------------------
/*!
  @brief Second-order correlation potential, Sigma, from Goldstone diagrams.

  @details
  Forms the second-order correlation potential for a valence state, as a
  non-local matrix on a radial sub-grid (every stride-th point of the full
  grid), by direct summation over a finite basis of internal core (hole) and
  excited states. Direct and exchange diagrams, with optional effective
  screening (fk) and hole-particle (etak) factors, and the two-body Breit
  correction.

  Used by CorrelationPotential: for all diagrams (Goldstone method), or for
  the exchange and Breit diagrams only (Feynman method).
*/
class Goldstone {

  // Pointer to shared radial grid (full grid)
  std::shared_ptr<const Grid> m_grid;

  using Basis = std::vector<DiracSpinor>;
  std::pair<Basis, Basis> m_basis;
  Coulomb::YkTable m_Yeh;

  // Parameters of the sub-grid: initial/final points, stride
  std::size_t m_i0, m_stride, m_subgrid_points;

  int m_n_min_core;
  bool m_include_G;

  std::optional<HF::Breit> m_Br{};

public:
  /*!
    @brief Splits basis into internal holes and excited states; forms Yk table.

    @details
    Basis states found in @p core with n >= @p n_min_core are the internal
    holes; basis states not in @p core are the excited states. Calculates the
    Yk integrals for each excited-hole pair.

    @param basis       Basis for internal lines (core and excited states);
                       must not be empty
    @param core        Core states: used to split @p basis into holes and
                       excited states
    @param i0          Full-grid index of the first sub-grid point
    @param stride      Sub-grid stride: every stride-th point of the full grid
    @param size        Number of sub-grid points
    @param n_min_core  Minimum n of core states included as internal holes
                       (lower core states are frozen)
    @param include_G   Include small-component (fg, gf, gg) parts of Sigma
    @param Br          Breit operator, for dSigma_Breit2 (copied; nullptr:
                       none)
  */
  Goldstone(const std::vector<DiracSpinor> &basis,
            const std::vector<DiracSpinor> &core, std::size_t i0,
            std::size_t stride, std::size_t size, int n_min_core = 1,
            bool include_G = false, const HF::Breit *Br = nullptr);

  /*!
    @brief Calculate Direct part of correlation potential.

    @details
    Direct diagrams, summed over internal holes and excited states, at energy
    @p en_v. Each multipolarity k is scaled by fk[k] * etak[k].

    @param kappa_v     Dirac quantum number of the valence state
    @param en_v        Energy at which Sigma is formed
    @param fk          Effective screening factors, for each multipole k
                       (missing: 1)
    @param etak        Effective hole-particle factors, for each k (missing:
                       1)
    @param n_max_core  Maximum n of internal hole states (<= 0: all)

    @return Direct part of Sigma, on the sub-grid.
  */
  GMatrix Sigma_direct(int kappa_v, double en_v,
                       const std::vector<double> &fk = {},
                       const std::vector<double> &etak = {},
                       int n_max_core = 99) const;

  /*!
    @brief Calculate Exchange part of correlation potential.

    @details
    @p fk are effective screening factors, applied to the multipole k of the
    outer Coulomb line; with @p fk_both_lines, also to the multipole of the
    inner line (inside \f$ P^k \f$), i.e., \f$ f_k f_l \f$.

    @return Exchange part of Sigma, on the sub-grid.
  */
  GMatrix Sigma_exchange(int kappa_v, double en_v,
                         const std::vector<double> &fk = {},
                         bool fk_both_lines = false) const;

  /*!
    @brief Calculate both parts (direct + exchange) of correlation potential.

    @details
    Direct and exchange diagrams formed together, in a single sum. Parameters
    as for Sigma_direct() and Sigma_exchange(); here, @p n_max_core applies to
    both the direct and exchange diagrams.

    @return Direct plus exchange parts of Sigma, on the sub-grid.
  */
  GMatrix Sigma_both(int kappa_v, double en_v,
                     const std::vector<double> &fk = {},
                     const std::vector<double> &etak = {}, int n_max_core = 99,
                     bool fk_both_lines = false) const;

  /*!
    @brief Calculates 2-body Breit correction to correlation potential.

    @details
    Second-order diagrams with one Coulomb line replaced by the two-body
    Breit interaction (factor of 2: either line), including direct and
    exchange parts. Each multipolarity k is scaled by fk[k]; etak[k] scales
    the direct part only.

    @param kappa_v        Dirac quantum number of the valence state
    @param en_v           Energy at which Sigma is formed
    @param fk             Effective screening factors, for each multipole k
    @param etak           Effective hole-particle factors, for each k (direct
                          part only)
    @param n_max_core     Maximum n of internal hole states (<= 0: all)
    @param max_n_breit    Maximum n of internal excited states (<= 0: all)

    @return Breit correction to Sigma, on the sub-grid.

    @note Must have Breit: returns zero if no Breit operator was given to the
          constructor.
  */
  GMatrix dSigma_Breit2(int kappa_v, double en_v,
                        const std::vector<double> &fk = {},
                        const std::vector<double> &etak = {},
                        int n_max_core = 99, int max_n_breit = -1) const;

  //! Internal {holes, excited} states (holes: core states with n >= n_min)
  const std::pair<Basis, Basis> &basis() const { return m_basis; }
  //! Yk table for each excited-hole pair
  const Coulomb::YkTable &Yeh() const { return m_Yeh; }

  //! Sub-grid stride: every stride-th point of the full grid
  std::size_t stride() const { return m_stride; }
  //! Minimum n of core states included as internal holes
  int n_min() const { return m_n_min_core; }
  //! Maximum l of the internal excited states
  int lmax() const { return DiracSpinor::max_l(m_basis.second); }

private:
  inline double get_k(int k, const std::vector<double> &f) const {
    const auto sk = std::size_t(k);
    return sk < f.size() ? f[sk] : 1.0;
  }
};

} // namespace MBPT