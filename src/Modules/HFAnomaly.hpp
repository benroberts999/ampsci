#pragma once
#include <array>
#include <string>
#include <utility>
class Wavefunction;
class DiracSpinor;
namespace IO {
class InputBlock;
}

namespace Module {

//! Calculates hyperfine anomaly: BW effect, differential anomaly, fits Rmag
void HFAnomaly(const IO::InputBlock &input, const Wavefunction &wf);

void b_plot(const IO::InputBlock &input, const Wavefunction &wf);

std::pair<std::array<double, 3>, std::array<double, 3>>
b_moments(const std::string &iso, const DiracSpinor &v, double R0_fm,
          int max_power);

} // namespace Module
