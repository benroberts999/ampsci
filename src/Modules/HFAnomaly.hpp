#pragma once
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Calculates hyperfine anomaly: BW effect, differential anomaly, fits Rmag
void HFAnomaly(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
