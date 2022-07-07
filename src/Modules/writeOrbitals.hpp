#pragma once
#include <string>
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Module: writes orbitals to text file (gnuplot format)
void writeOrbitals(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module