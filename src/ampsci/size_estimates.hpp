#pragma once
#include "IO/InputBlock.hpp"
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

/*!
  @brief Sizes of the main tables and matrices a run would form, estimated
  from its input file without running it.
  @details
  Counts are exact for the given input (they use the same selection rules and
  defaults as the run); the memory is an estimate (hash-table overhead for the
  integral tables; workspace for the dense CI diagonalisation).
  See @ref estimate_sizes.
*/
struct SizeEstimates {
  //! Size of one CI (J, parity) block
  struct CIBlock {
    //! Twice the total angular momentum, 2J
    int twoJ;
    //! Parity, +1 or -1
    int parity;
    //! Number of CSFs
    std::size_t num_CSFs;
    //! Peak memory while forming and diagonalising the dense CI matrix (GB)
    double matrix_GB;
    //! Size of this block in the ci solutions file (GB)
    double file_GB;
  };

  //! Full single-particle basis (Basis{} states) and core
  std::string basis_string{};
  std::string core_string{};
  std::size_t num_orbitals{0};
  std::size_t num_core{0};

  //! CI space
  std::string ci_basis_string{};
  std::size_t num_ci_orbitals{0};
  int num_electrons{2};

  //! All Coulomb Q^k integrals over the full basis, all k (as needed by
  //! structure radiation)
  std::size_t qk_all{0};
  //! Coulomb Q^k integrals among the CI orbitals
  std::size_t qk_CI{0};
  //! Coulomb Q^k integrals with one or two core orbitals (for MBPT)
  std::size_t qk_MBPT{0};
  //! Whether this run forms each of the Q^k tables
  bool qk_all_required{false};
  bool qk_CI_required{false};
  bool qk_MBPT_required{false};
  //! Breit B^k integrals (zero unless this run forms them)
  std::size_t bk{0};

  //! One entry per requested (J, parity)
  std::vector<CIBlock> ci_blocks{};
};

/*!
  @brief Estimates the sizes of the integral tables and CI matrices of a run,
  from its input.
  @details
  Uses the HartreeFock{} core, the Basis{} states, and the CI{} options (with
  the same defaults as the CI run). Orbitals are labels only: no
  wavefunctions are calculated. Counted:
  - all Q^k over the full basis, for all k (the table structure radiation
    uses; required if structure radiation is requested);
  - Q^k among the CI orbitals (required if there is a CI{} block);
  - Q^k with one or two core orbitals, over the MBPT basis (required when
    Sigma_1 is calculated in the CI - no Correlations{} block, or a
    dSigma/dE correction that the correlation potential does not supply -
    or Sigma_2 is included);
  - B^k over Breit_basis (if Breit is included in HF);
  - the number of CSFs, peak memory, and ci-file size of each requested
    (J, parity), from @ref CI::requested_J_pi.

  @param input  The full ampsci input.
  @return The estimates. Without a CI{} block, the CI and MBPT counts use the
          CI defaults, and are not marked required.

  @note Counting the all-Q^k and MBPT integrals over a large basis takes
        noticeable time (it loops over all sets of four orbitals).
*/
[[nodiscard]] SizeEstimates estimate_sizes(const IO::InputBlock &input);

/*!
  @brief Prints the estimates of @ref estimate_sizes.
  @details
  All tables are printed; those the run requires are marked, and counted in
  the total, with the largest CI block (blocks are solved one at a time).
*/
void print_size_estimates(const SizeEstimates &estimates,
                          std::ostream &os = std::cout);
