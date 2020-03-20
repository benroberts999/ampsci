#pragma once
#include <memory>
#include <string>
#include <vector>
class Wavefunction;
class UserInputBlock;
namespace DiracOperator {
class TensorOperator;
}

namespace Module {

//! Calculates matrix elements of any tensor operator, with RPA
void matrixElements(const UserInputBlock &input, const Wavefunction &wf);

void calculateBohrWeisskopf(const UserInputBlock &input,
                            const Wavefunction &wf);

//! Calculates state lifetimes (using E1 and E2 only). nb: HF energies
void calculateLifetimes(const UserInputBlock &input, const Wavefunction &wf);

//! Module: Calculates 2nd order energy shift for valence states
void SecondOrder(const UserInputBlock &input, const Wavefunction &wf);

//! Returns a ptr to the requested operator, with given properties
std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(const std::string &operator_str, const UserInputBlock &input,
                 const Wavefunction &wf);

} // namespace Module
