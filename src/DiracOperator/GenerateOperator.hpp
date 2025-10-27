#pragma once
#include "IO/InputBlock.hpp"
#include "Operators.hpp"
#include "TensorOperator.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <memory>
#include <string>
#include <vector>

namespace DiracOperator {

//! List of avilable operators, and their generator functions
static const std::vector<std::pair<
    std::string, std::unique_ptr<DiracOperator::TensorOperator> (*)(
                     const IO::InputBlock &input, const Wavefunction &wf)>>
    operator_list{{"E1", &generate_E1},
                  {"E1v", &generate_E1v},
                  {"E2", &generate_E2},
                  {"ialpha", &generate_ialpha},
                  {"Ek", &generate_Ek},
                  {"Ek_omega", &generate_Ek_omega},
                  {"Mk_omega", &generate_Mk_omega},
                  {"Ekv_omega", &generate_Ekv_omega},
                  {"M1", &generate_M1},
                  {"M1nr", &generate_M1nr},
                  {"hfs", &generate_hfs},
                  {"fieldshift", &generate_fieldshift},
                  {"r", &generate_r},
                  {"sigma_r", &generate_sigma_r},
                  {"pnc", &generate_pnc},
                  {"Vrad", &generate_Vrad},
                  {"MLVP", &generate_MLVP},
                  {"dr", &generate_dr},
                  {"p", &generate_p},
                  {"l", &generate_l},
                  {"s", &generate_s}};

//------------------------------------------------------------------------------

//! Returns a unique_ptr (polymorphic) to the requested operator, with given
//! properties
std::unique_ptr<DiracOperator::TensorOperator>
generate(std::string_view operator_name, const IO::InputBlock &input,
         const Wavefunction &wf);

//! List available operators
void list_operators();

} // namespace DiracOperator
