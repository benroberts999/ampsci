#pragma once
#include <iostream>
class Wavefunction;
namespace IO {
class InputBlock;
}

namespace Module {

//! Calculates atomic ionisation factor, Kion(E,q) = |<a|j_L|e>|^2, summed over
//! core states, where E is deposited energy: E = |E_a| + e
void Kionisation(const IO::InputBlock &input, const Wavefunction &wf);

void photo(const IO::InputBlock &input, const Wavefunction &wf);

void formFactors(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
