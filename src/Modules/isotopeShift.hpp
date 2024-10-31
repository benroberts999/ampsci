#pragma once

// Forward declare classes:
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

void fieldShift(const IO::InputBlock &input, const Wavefunction &wf);

void fieldShift_direct(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module