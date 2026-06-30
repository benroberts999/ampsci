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
  Sigma_L, for each valence state, and writes these to the .sl file.
  The .sl file may then be read in by the Correlations block (ladder_file
  option) to include Sigma_L into the correlation potential.
  Runs before Correlations. Options are parsed from the Ladder{} input block.
*/
void ladder(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace MBPT
