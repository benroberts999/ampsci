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
  ExternalField::TDHF S(&d, wf.vHF());
  ExternalField::TDHF T(&d, wf.vHF());
  ExternalField::TDHF_DCP ddcp(&S, &T, 2, wf.vHF());
  std::cout << __LINE__ << std::endl;

  S.solve_core(0.0);

  T.solve_core(0.0);

  ddcp.solve_core(0.0);

  //
}

} // namespace Module
