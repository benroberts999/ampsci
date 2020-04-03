#pragma once
class Wavefunction;
namespace IO {
class UserInputBlock;
class UserInput;
} // namespace IO

namespace Module {

//! Calculates PNC using sum-over-states and solving-equations method (HF+RPA)
void calculatePNC(const IO::UserInputBlock &input, const Wavefunction &wf);

void polarisability(const IO::UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
