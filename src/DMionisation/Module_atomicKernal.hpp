#pragma once
#include <iostream>
class Wavefunction;
class UserInputBlock;

namespace Module {

void atomicKernal(const UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
