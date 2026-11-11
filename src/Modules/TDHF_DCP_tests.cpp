#include "DiracOperator/include.hpp" //For E1 operator
#include "ExternalField/TDHF_DCP.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

// Declare, register, then define below.
void TDHF_DCP_tests(const IO::InputBlock &input, const Wavefunction &wf);
namespace {
const Register r_TDHF_DCP_tests{
  "TDHF_DCP_tests", "Module for testing DCP by TDHF", &TDHF_DCP_tests};
} // namespace

void TDHF_DCP_tests(const IO::InputBlock &input, const Wavefunction &wf) {
  (void)input; // suppress unused variable warning

  const auto d = DiracOperator::E1(wf.grid());
  ExternalField::TDHF S(&d, wf.vHF());
  ExternalField::TDHF T(&d, wf.vHF());
  ExternalField::TDHF_DCP ddcp(&T, &S, 2, wf.vHF());
  std::cout << __LINE__ << std::endl;

  S.solve_core(0.0);

  T.solve_core(0.0);

  ddcp.solve_core(0.0);

  //
}

} // namespace Module
