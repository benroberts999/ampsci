#pragma once
#include "IO/InputBlock.hpp"
#include "TensorOperator.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <memory>
#include <string>
#include <vector>

namespace DiracOperator {

//------------------------------------------------------------------------------

//! Returns a unique_ptr (polymorphic) to the requested operator, with given
//! properties
std::unique_ptr<DiracOperator::TensorOperator>
generate(std::string_view operator_name, const IO::InputBlock &input,
         const Wavefunction &wf);

//! List available operators
void list_operators();

} // namespace DiracOperator
