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

void vertexQED(const IO::InputBlock &input, const Wavefunction &wf);

//! Used for finding A and b for effective vertex QED operator
void hyperfine_vertex_test(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
