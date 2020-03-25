#pragma once
class Wavefunction;
class UserInputBlock;
class UserInput;

namespace Module {

//! Calculates PNC using sum-over-states and solving-equations method (HF+RPA)
void calculatePNC(const UserInputBlock &input, const Wavefunction &wf);

} // namespace Module
