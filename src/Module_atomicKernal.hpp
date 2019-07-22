#pragma once
#include <iostream>
class Wavefunction;
class UserInput;

namespace Module {

static const std::string ThisModule = "Module::AtomicKernal";

void atomicKernal(const UserInput &input, const Wavefunction &wf);

} // namespace Module
