#include "DiracOperator/include.hpp" //For E1 operator
#include "ExternalField/TDHF_DCP.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/exampleModule.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void TDHF_DCP_tests(const IO::InputBlock &input, const Wavefunction &wf) {
  (void)input; // suppress unused variable warning

  const auto d = DiracOperator::E1(wf.grid());

  ExternalField::TDHF_DCP ddcp(&d, &d, 2, wf.vHF());

  ddcp.solve_core(0.0);

  //
}

} // namespace Module
