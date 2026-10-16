#pragma once
#include "CI_Integrals.hpp"
#include "CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>
#include <string>
#include <vector>

/*! 
  @brief Functions and classes for Configuration Interaction calculations

  @details
    Main functions are: 
    - @ref configuration_interaction
    - @ref run_CI
    
    Main Classes are: 
    - @ref CSF2, @ref CSF3
    - @ref PsiJPi
*/
namespace CI {

/*!
  @brief Runs Configuration Interation: returns CI solutions for all
  requested J and parity values.

  @details
  Reads options from @p input, the CI Input Block, and constructs the CI basis
  from @p wf, computes the required Coulomb (and optionally Breit and two-body 
  MBPT) integrals, then calls run_CI() for each requested (J, parity) pair.

  The returned vector contains one @ref PsiJPi per {J, parity} combination, each
  holding the eigenvalues and CI expansion coefficients for the requested number
  of solutions.

  }
  @endcode
  * Always check for up-to-date options from command line: `$ ampsci -i CI`
  * See also @ref run_CI, which this function calls

  @param input   Input block containing CI options.
  @param wf      Fully initialised Wavefunction object supplying the orbital
                 basis and radial grid.
  @return @ref Solutions: one PsiJPi per (J, parity) combination requested,
          together with the integral tables used to build the CI Hamiltonians.

  @note If run with `read_only`, no integrals are calculated, and the returned
        integral tables are empty (see Integrals::availableQ()).
*/
Solutions configuration_interaction(const IO::InputBlock &input,
                                    const Wavefunction &wf);

/*!
  @brief The (2J, parity) blocks requested in a CI input block.
  @details
  'twoJ', 'twoJ+', 'twoJ-' (lists of integer 2J) take precedence over 'J',
  'J+', 'J-' (integers, fractions such as 3/2, or decimals such as 1.5). A
  parity-specific list takes precedence over the plain one, for its parity.
  2J must be even for two electrons and odd for three; other values are
  skipped, with a warning (once per value).

  @param input          The CI input block.
  @param num_electrons  Number of valence electrons (2 or 3).
  @return List of {2J, parity}, even parity first.
*/
std::vector<std::pair<int, int>> requested_J_pi(const IO::InputBlock &input,
                                                int num_electrons);

/*!
  @brief Constructs and solves the CI eigenvalue problem for a single J,pi
  @details
  Builds the CI+MBPT Hamiltonian matrix in the basis of two- or three-electron
  configuration state functions (CSFs; see @ref CSF2, @ref CSF3) with total
  angular momentum @p twoJ /2 and parity @p parity, 
  then solves the eigenvalue problem to obtain CI energies and
  expansion coefficients.

  The Hamiltonian includes:
  - one-body terms from @p h1
  - two-body Coulomb interaction via the \f$ Q^k \f$ integrals in @p qk
  - optionally, two-body Breit corrections via \f$ B^k \f$ integrals in @p Bk
    (used if @p Bk is non-empty)
  - optionally, two-body MBPT \f$ \Sigma_2 \f$ corrections via \f$ S^k \f$
    integrals in @p Sk (used if @p include_Sigma2 is true, and @p Sk is non-empty)

  The number of solutions returned is controlled by @p num_solutions and
  @p all_below: if @p all_below is set it takes precedence and all eigenstates
  with total energy below the threshold are found.

  @param ci_sp_basis  Single-particle basis states spanning the CI space.
  @param twoJ         Twice the total angular momentum, 2J (non-negative; even
                      for two electrons, odd for three).
  @param parity       Parity of the block: +1 (even) or -1 (odd).
  @param num_electrons Number of valence electrons: 2 or 3.
  @param num_solutions Number of lowest eigenstates to find. Ignored if
                       @p all_below_cm is set. Pass 0 to find all solutions.
  @param all_below_cm If set, find all eigenstates with total energy below this
                      value (in cm^-1). Overrides @p num_solutions.
  @param h1           Table of one-body Hamiltonian matrix elements between
                      single-particle basis states.
  @param qk           Table of two-body Coulomb \f$ Q^k \f$ integrals.
  @param Bk           Table of two-body Breit \f$ B^k \f$ integrals. Ignored
                      (treated as absent) if the table is empty.
  @param Sk           Table of two-body MBPT \f$ \Sigma_2 \f$ (\f$ S^k \f$)
                      integrals. Only used when @p include_Sigma2 is true.
  @param include_Sigma2 If true, add two-body MBPT corrections from @p Sk to
                        the CI Hamiltonian.
  @param print_details  If true, print a breakdown of the leading configurations
                        for each solution. Leads to very large output if
                        @p num_solutions is large
  @param read_only    If true, only the solutions already in the file are read;
                      the eigenvalue problem is not solved, and nothing is
                      written [default: false].
  @param outstream    Output stream for progress and results [default: stdout].
  @param ci_fname     Filename for reading/writing CI solutions ("" disables).
  @param s1c          Pointer to derivative (dSigma/dE) correction for
                      Sigma_1; ignored if nullptr. See @ref Sigma1Correction.
                      Its reference energy E0 is the lowest energy of this
                      (J, parity), so is found self-consistently: see
                      @ref iterate_E0.
  @param hk           Average S^k/Q^k ratios, for extrapolating Sigma_2 beyond
                      the cis2 basis; empty for no extrapolation.
                      See @ref MBPT::average_hk.
  @param dSk          Pointer to the dS^k/dE0 table (Brillouin-Wigner Sigma_2);
                      ignored if nullptr. See @ref CI::corrected_Sk.
  @param E0_sigma2    Reference E0 that Sk and dSk were tabulated at.
                      
  @return PsiJPi (@ref PsiJPi) containing the CI eigenvalues and expansion
  coefficients for the requested solutions.
*/
PsiJPi run_CI(const std::vector<DiracSpinor> &ci_sp_basis, int twoJ, int parity,
              int num_electrons, int num_solutions,
              std::optional<double> all_below_cm,
              const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
              const Coulomb::WkTable &Bk, const Coulomb::LkTable &Sk,
              bool include_Sigma2, bool print_details, bool read_only = false,
              std::ostream &outstream = std::cout,
              const std::string &ci_fname = "",
              const Sigma1Correction *s1c = nullptr,
              const std::vector<double> &hk = {},
              const Coulomb::LkTable *dSk = nullptr, double E0_sigma2 = 0.0);

} // namespace CI