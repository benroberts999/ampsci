#pragma once
#include "DiracOperator/include.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "Modules/exampleModule.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void TDHF_DCP_tests(const IO::InputBlock &input, const Wavefunction &wf);

} // namespace Module
