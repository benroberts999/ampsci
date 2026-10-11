#pragma once

// Forward declarations
namespace IO {
class InputBlock;
}
class Wavefunction;

namespace MBPT {

/*!
  @brief Driver for the ladder-diagram calculation: the Ladder{} input block.

  @details
  Calculates (and iterates to convergence) the ladder integrals Lk, writing
  them to the .lk file. Then constructs the ladder correlation potential,
  Sigma_L, for the lowest valence state of each kappa (or for each valence
  state, if each_valence is set), and writes these to the .sl file.
  The .sl file may then be read in by the Correlations block (ladder
  option) to include Sigma_L into the correlation potential.
  Runs before Correlations. Options are parsed from the Ladder{} input block.

  @param input  The Ladder{} input block
  @param wf     Wavefunction: supplies the core, valence states, and basis

  @warning: Experimental, not fully tested, be careful

  @note Valence states missing from the basis are skipped (with a warning).
*/
void ladder(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace MBPT
