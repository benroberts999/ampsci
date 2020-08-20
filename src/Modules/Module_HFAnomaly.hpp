#pragma once
#include <memory>
#include <string>
#include <vector>
class Wavefunction;
namespace IO {
class UserInputBlock;
}
namespace DiracOperator {
class TensorOperator;
}

namespace Module {

//! Calculates Bohr-Weisskopf effect for hyperfine structure
void calculateBohrWeisskopf(const IO::UserInputBlock &input,
                            const Wavefunction &wf);

//! Calculates hyperfine anomaly for list of isotopes
void HFAnomaly(const IO::UserInputBlock &input, const Wavefunction &wf);

void HF_rmag(const IO::UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
