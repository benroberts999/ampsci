#pragma once
#include <iostream>
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

void atomicKernel(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
