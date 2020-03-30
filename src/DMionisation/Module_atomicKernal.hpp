#pragma once
#include <iostream>
class Wavefunction;
namespace IO {
class UserInputBlock;
}

namespace Module {

void atomicKernal(const IO::UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
