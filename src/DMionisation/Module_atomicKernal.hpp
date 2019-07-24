#pragma once
#include <iostream>
class Wavefunction;
class UserInputBlock;

namespace Module {

static const std::string ThisModule = "Module::AtomicKernal";

void atomicKernal(const UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
