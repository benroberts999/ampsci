#pragma once
#include "IO/InputBlock.hpp"
#include "Operators/include.hpp"
#include "TensorOperator.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <memory>
#include <string>
#include <vector>

namespace DiracOperator {

/*!
  @brief List of available operators, their generator functions, and a short description.
  @details
  This list is used by @ref generate() to form operators, and is called in a few modules.

  It is also what gets printed to the screen by:
  ```shell
    ampsci -o
  ```

  @note Operators must be added to this list in order for ampsci to know where to find them, and to document to users what's available.
*/
static const std::vector<
  std::tuple<std::string,
             std::unique_ptr<DiracOperator::TensorOperator> (*)(
               const IO::InputBlock &input, const Wavefunction &wf),
             std::string>>
  operator_list{
    {"E1", &generate_E1, "Electric dipole (moment), length form: -|e|r"},
    {"E1v", &generate_E1v, "Electric dipole, v-form"},
    {"E2", &generate_E2, "Electric quadrupole moment operator"},
    {"Ek", &generate_Ek, "Electric multipole moment operator, in low qr limit"},
    {"ialpha", &generate_ialpha, "i*alpha (propto E1v)"},
    {"M1", &generate_M1, "Magnetic dipole (relativistic formula)"},
    {"M1nr", &generate_M1nr, "Non-relativistic M1"},
    {"Multipole", &generate_Multipole,
     "Multipole transition operators (Vector,Axial,Scalar,Pseudoscalar)"},
    {"hfs", &generate_hfs, "Hyperfine structure k-pole operators"},
    {"fieldshift", &generate_fieldshift, "Field-shift F(r) operator"},
    {"r", &generate_r, "radial (scalar) |r|"},
    {"sigma_r", &generate_sigma_r, "scalar sigma.r operator"},
    {"pnc", &generate_pnc, "NSI PNC operator"},
    {"Vrad", &generate_Vrad, "QED Radiative potential"},
    {"MLVP", &generate_MLVP,
     "Magnetic-Loop vacuum polarisation vertex correction to HFS"},
    {"p", &generate_p, "Momentum operator"},
    {"l", &generate_l, "Orbital L"},
    {"s", &generate_s, "Spin S (not sigma)"}};

//------------------------------------------------------------------------------

//! Returns a unique_ptr (polymorphic) to the requested operator, with given properties
/*!
  @details

  From the command line, use
  ```shell
    ampsci -o
  ```
  to get a list of available operators.
  For a specific operator 'OperatorName', use:
  ```shell
    ampsci -o OperatorName
  ```
  to see available run-time options for that operator.

*/
std::unique_ptr<DiracOperator::TensorOperator>
generate(std::string_view operator_name, const IO::InputBlock &input,
         const Wavefunction &wf);

//! List available operators
void list_operators();

} // namespace DiracOperator
