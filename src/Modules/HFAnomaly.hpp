#pragma once
#include <memory>
#include <string>
#include <vector>
class Wavefunction;
namespace IO {
class InputBlock;
}
namespace DiracOperator {
class TensorOperator;
}

namespace Module {

//! Calculates Bohr-Weisskopf effect for hyperfine structure, including screening factors
void BohrWeisskopf(const IO::InputBlock &input, const Wavefunction &wf);

//! Calculates hyperfine anomaly for list of isotopes
/*! @details
```cpp
Module::HFAnomaly{
  rpa; // true/false include RPA? (uses diagram method)
  options; //options for HFS operator (typically left blank)
  A; // List of A's, comma-separated (uses wavefunction A as reference)
}
```
*/
void HFAnomaly(const IO::InputBlock &input, const Wavefunction &wf);

//! Finds magnetic radius that reproduced given hyperfine anomaly
/*! @details
  For isotope 1 and 2
  Loops over many values for magnetic radius of isotope 1, Rmag(1).
  For each, finds Rmag(2) that reproduces a given hyperfine anomaly.
  Potential to use two states to find rmag. If not, still interesting?
*/
void HF_rmag(const IO::InputBlock &input, const Wavefunction &wf);

//! Calculates eta_sp, which links BW effect is s tp p states of H-like ions
void BW_eta_sp(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
