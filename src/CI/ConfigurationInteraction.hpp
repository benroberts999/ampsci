#pragma once
#include "CSF.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/meTable.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>
#include <vector>

//! Functions and classes for Configuration Interaction calculations
/*! @details
    Main functions are: 
    - @ref configuration_interaction
    - @ref run_CI
    
    Main Classes are: 
    - @ref CSF2
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

  As of writing, options are:
  @code{.java}
  // Available CI options/blocks
  CI{
    ci_basis;
      // Basis used for CI expansion; must be a sub-set of
      // full ampsci basis [default: 20spdf]
    J;
      // List of total angular momentum J for CI solutions
      // (comma separated). Must be integers (two-electron
      // only). []
    J+;
      // As above, but for EVEN CSFs only (takes precedence
      // over J).
    J-;
      // As above, but for ODD CSFs (takes precedence over J).
    num_solutions;
      // Number of CI solutions to find (for each J/pi) [5]
    all_below_cm;
      // Find all CI solutions for energies below this
      // threshold, in inverse cm. Note that this is the total
      // energy, not the excitation energy. If set,
      // num_solutions is ignored.
    sigma1;
      // Include one-body MBPT correlations? [false]
    sigma2;
      // Include two-body MBPT correlations? [false]
    Brueckner;
      // Use Brueckner (spectrum) states for CI basis? Must
      // have Spectrum and sigma1. [false]
    cis2_basis;
      // The subset of ci_basis for which the two-body MBPT
      // corrections are calculated. Must be a subset of
      // ci_basis. If existing sk file has more integrals,
      // they will be used. [default: Nspdf, where N is
      // maximum n for core + 3]
    Breit2;
      // Include two-body Breit? Default is true if Breit
      // included in HF. Ignored if Breit not included in HF.
      // [true]
    Breit_basis;
      // Subset of ci_basis used to include two-body Breit
      // corrections into CI matrix. Large basis is slow, uses
      // huge memory, and makes small contribution. [default:
      // Nspdf, where N is maximum n for core + 6]
    s1_basis;
      // Usually should be left as default. Basis used for the
      // one-body MBPT diagrams (Sigma^1) internal lines.
      // These are the most important, so in general the
      // default (all basis states) should be used. Must be a
      // subset of full ampsci basis. [default: full basis]
      //  - Note: if CorrelationPotential is available, it
      // will be used instead of calculating the Sigma_1
      // integrals
    s2_basis;
      // Usually should be left blank. Basis used for internal
      // lines of the two-body MBPT diagrams (Sigma^2)
      // internal lines. Must be a subset of s1_basis.
      // [default: s1_basis]
    n_min_core;
      // Minimum n for core to be included in MBPT [1]
    max_k;
      // Maximum k (multipolarity) to include when calculating
      // new Coulomb integrals. Higher k often contribute
      // negligably. Note: if qk file already has higher-k
      // terms, they will be included. Set negative (or very
      // large) to include all k. [8]
    denominators;
      // 'RS', 'Fermi', 'Fermi0'. Denominators used in Sigma2
      // matrix elements. RS uses actual states for external
      // legs, Fermi uses the lowest excited state for each
      // kappa, Fermi0 uses lowest excited state for all
      // kappas (and thus cancels in all except diagram 'd').
      // [Fermi0]
    qk_file;
      // Filename for storing two-body Coulomb integrals. By
      // default, is ~ At.qk, where At is atomic symbol +
      // 'identity'.
    sk_file;
      // Filename for storing two-body Sigma_2 integrals. By
      // default, is At_n_b_k.sk, where At is atomic symbol, n
      // is n_min_core, b is cis2_basis, k is max_k.
    bk_file;
      // Filename for storing two-body Breit integrals. By
      // default, is ~ At.bk, where At is atomic symbol +
      // 'identity'.
    no_new_integrals;
      // Usually false. If set to true, ampsci will not
      // calculate any new Coulomb or Sigma_2 integrals, even
      // if they are implied by the above settings. This saves
      // time when we know all required integrals already
      // exist, since the code doesn't need to check. [true]
    exclude_wrong_parity_box;
      // Excludes the Sigma_2 box corrections that have
      // 'wrong' parity when calculating Sigma2 matrix
      // elements. Note: If existing sk file already has
      // these, they will be included [false]
    sort_output;
      // Sort output by energy? Default is to sort by J and Pi
      // first. [false]
    print_details;
      // Condition to print details of each CI solution
      // (otherwise just prints summary) [true]
    parallel_ci;
      // Run CI in parallel (solve each J/Pi in parallel).
      // Faster, uses slightly more memory [true]
  }
  @endcode
  * Always check for up-to-date options from command line: `$ ampsci -i CI`
  * See also @ref run_CI, which this function calls

  @param input   Input block containing CI options.
  @param wf      Fully initialised Wavefunction object supplying the orbital
                 basis and radial grid.
  @return Vector of PsiJPi, one entry per (J, parity) combination requested.
*/
std::vector<PsiJPi> configuration_interaction(const IO::InputBlock &input,
                                              const Wavefunction &wf);

/*!
  @brief Constructs and solves the CI eigenvalue problem for a single J,pi
  @details
  Builds the CI+MBPT Hamiltonian matrix in the basis of two-electron 
  configuration state functions (CSFs) with total angular momentum 
  @p twoJ /2 and parity @p parity, 
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
  @param twoJ         Twice the total angular momentum, 2J (must be a
                      non-negative even integer for two-electron systems).
  @param parity       Parity of the sector: +1 (even) or -1 (odd).
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
  @param outstream    Output stream for progress and results [default: stdout].
  @return PsiJPi (@ref PsiJPi) containing the CI eigenvalues and expansion 
  coefficients for the requested solutions.
*/
PsiJPi run_CI(const std::vector<DiracSpinor> &ci_sp_basis, int twoJ, int parity,
              int num_solutions, std::optional<double> all_below_cm,
              const Coulomb::meTable<double> &h1, const Coulomb::QkTable &qk,
              const Coulomb::WkTable &Bk, const Coulomb::LkTable &Sk,
              bool include_Sigma2, bool print_details,
              std::ostream &outstream = std::cout);

} // namespace CI