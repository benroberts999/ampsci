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

//! Calculates Bohr-Weisskopf effect for hyperfine structure
void calculateBohrWeisskopf(const IO::InputBlock &input,
                            const Wavefunction &wf);

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

void HF_rmag(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
