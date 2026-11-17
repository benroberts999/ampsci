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
  With include_P, the ladder equation also includes P (screening, with its
  exchange and hole-particle parts); only the valence energy corrections are
  printed, without the part the Feynman direct diagram contains, and Sigma_L
  is not formed. With include_R (needs include_P), the ladder equation also
  includes the core relaxation R (the core single excitations rho_ma), and
  the interaction with the change of the core density, <v|dV^R|v>, is added
  to the valence energy corrections.

  @param input  The Ladder{} input block
  @param wf     Wavefunction: supplies the core, valence states, and basis

  @warning: Experimental, not fully tested, be careful

  @note Valence states missing from the basis are skipped (with a warning).
*/
void ladder(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace MBPT
