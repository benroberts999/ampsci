#pragma once
#include "Angular/SixJTable.hpp"
#include "Coulomb/QkTable.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <optional>
#include <string>
#include <vector>

// Forward declarations
namespace IO {
class InputBlock;
}
class Wavefunction;

namespace MBPT {

/*!
  @brief Options controlling the ladder-diagram correlation potential.
  @details
  Parsed from the Correlations{ Ladder{} } sub-block. The ladder is built on its
  own basis subset and n_min_core (independent of the rest of Sigma); only the
  sub-grid is shared with the correlation potential. Qk and Lk are cached to disk
  and read back when present.

  @note basis_str selects the excited-state subset (CI::basis_subset). Empty =>
        entire excited basis.
*/
struct LadderOptions {
  std::string basis_str{};
  int n_min_core{1};
  int max_k{8};
  bool include_L4{false};
  bool full_basis{false};
  int max_it{15};
  double damp{0.0};
  double eps_target{1.0e-5};
  std::string Qk_file{};
  std::string Lk_file{};
  bool from_scratch{false};
};

/*!
  @brief Parses (and documents) the Correlations{ Ladder{} } sub-block.
  @details
  Reads the ladder options from the Ladder sub-block of @p correlations, filling
  default Qk/Lk cache filenames from @p identity. Returns std::nullopt if no
  Ladder sub-block is present (ladder disabled) or if help was requested.

  So that `ampsci -i Correlations` documents them, a 'help' option is propagated
  into the ladder block when @p correlations is in help mode (even if no Ladder
  sub-block is present).

  @param correlations  The Correlations input block
  @param identity       Atom identity, for default cache filenames
  @return Populated LadderOptions, or std::nullopt if ladder is off
*/
std::optional<LadderOptions>
parse_ladder_options(const IO::InputBlock &correlations,
                     const std::string &identity);

//------------------------------------------------------------------------------
/*!
  @brief Ladder-diagram correlation potential, Sigma_L.
  @details
  Companion to MBPT::Goldstone and MBPT::Feynman. Holds the Coulomb table Qk and
  the converged ladder table Lk (plus the hole/excited split and a 6j table), and
  builds the ladder correction to the correlation potential, Sigma_L(en), on a
  given sub-grid.

  The constructor reproduces the standalone ladder calculation: it builds (or
  reads) Qk, fills and iterates Lk to convergence, prints the MBPT(2)/ladder
  energy diagnostics, and (if full_basis) extends Qk to the full basis. The only
  difference from the standalone module is that the Sigma_L matrix is built later
  (see Sigma_ladder / CorrelationPotential::formSigma_L).

  @note Sub-grid (r0, rmax, stride) is supplied by the owning CorrelationPotential
        so Sigma_L lands on exactly the same grid as the base Sigma.
*/
class LadderPotential {
  std::vector<DiracSpinor> m_holes{};
  std::vector<DiracSpinor> m_excited{};
  std::vector<DiracSpinor> m_valence{}; // basis-version valence (Lk i-orbitals)
  // Projection basis for Sigma_L: the full basis when full_basis is set
  // (otherwise unused -- projection is onto the single valence state).
  std::vector<DiracSpinor> m_proj_basis{};
  Coulomb::QkTable m_qk{};
  Coulomb::LkTable m_lk{};
  Angular::SixJTable m_sjt{};
  double m_r0;
  double m_rmax;
  std::size_t m_stride;
  bool m_include_L4;
  bool m_full_basis;

public:
  /*!
    @brief Constructs the ladder potential: builds/reads Qk and iterates Lk.
    @param wf            Solved wavefunction (basis, core, valence, etc.)
    @param r0,rmax,stride  Sub-grid (from CorrelationPotential)
    @param options       Ladder options (basis subset, n_min_core, caching, ...)
  */
  LadderPotential(const Wavefunction &wf, double r0, double rmax,
                  std::size_t stride, const LadderOptions &options);

  //! Ladder correction Sigma_L for the given kappa/energy, on the sub-grid.
  //! @details Projection (Sigma_ladder filters by kappa internally): the full
  //! basis if full_basis was set, otherwise the single valence state Fv.
  //! include_G controls whether the lower (g) component is included.
  GMatrix Sigma_ladder(int kappa, double ev, const DiracSpinor *Fv,
                       bool include_G) const;
};

} // namespace MBPT
