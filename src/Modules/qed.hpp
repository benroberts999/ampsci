#pragma once
#include <string>
#include <vector>
class Wavefunction;
namespace IO {
class InputBlock;
}
namespace DiracOperator {
class TensorOperator;
}

namespace Module {

//! Calculates QED corrections to energies and matrix elements
//!@details: Note: Should be checked against manual calc, since not 100%
//! included
void QED(const IO::InputBlock &input, const Wavefunction &wf);

//! Calculates vertex QED corrections matrix elements
void vertexQED(const IO::InputBlock &input, const Wavefunction &wf);

std::vector<std::string> calc_vertexQED(const IO::InputBlock &input,
                                        const Wavefunction &wf);

} // namespace Module
