#include "Modules/CI_Pol.hpp"
#include "CI/CI_Integrals.hpp"
#include "DiracOperator/include.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void CI_Pol(const IO::InputBlock &input, const Wavefunction &wf) {
  // This module computes matrix elements and energies for two valence systems

  //these quantities are then used to compute second order amplitudes using the sum over states (SOS) method

  // In this example, we will solve a new wavefunction, assuming a different
  // nuclear charge distribution, and see the difference in the energies and E1
  // matrix elements this produces.

  // Read in some optional input options: A, rms, and type
  // "checkBlock" is optional, but helps prevent errors on input options:
  input.check({{"J", "Angular momentum of wavefunction"},
               {"parity", "Parity of state"},
               {"K", "Rank of operator"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto J = input.get("J", 0);
  const auto parity = input.get("parity", 1);

  const auto wfV = wf.CIwf(J, parity);

  // get allowed states for matrix elements of oprator of rank K

  const auto K = input.get("K", 0);
  for (int J2 = J - K; J2 <= J + K; J2++) {
    if (J2 < 0)
      continue;
    if (J2 == 0 && J == 0)
      continue;
    const auto wfW = wf.CIwf(J2, -parity);
    CI::ReducedME(wfV, 0, wfW, 0, table, K, )
  }

  const auto J = input.get("J", 0);
  const auto parity = input.get("parity", 1);

  const auto wfv = wf.CIwf(J, parity);
}
} // namespace Module
