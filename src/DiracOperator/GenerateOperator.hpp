#pragma once
#include "IO/InputBlock.hpp"
#include "Operators/include.hpp"
#include "TensorOperator.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <memory>
#include <string>
#include <vector>

namespace DiracOperator {

//! List of avilable operators, and their generator functions
static const std::vector<
    std::tuple<std::string,
               std::unique_ptr<DiracOperator::TensorOperator> (*)(
                   const IO::InputBlock &input, const Wavefunction &wf),
               std::string>>
    operator_list{
        {"E1", &generate_E1, "Electric dipole (moment), length form: -|e|r"},
        {"E1v", &generate_E1v, "Electric dipole, v-form"},
        {"E2", &generate_E2, "Electric quadrupole moment operator"},
        {"Ek", &generate_Ek,
         "Electric multipole moment operator, in low qr limit"},
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
        {"dr", &generate_dr, "Radial scalar derivative"},
        {"p", &generate_p, "Momentum operator"},
        {"l", &generate_l, "Orbital L"},
        {"s", &generate_s, "Spin S (not sigma)"}};

//------------------------------------------------------------------------------

//! Returns a unique_ptr (polymorphic) to the requested operator, with given
//! properties
std::unique_ptr<DiracOperator::TensorOperator>
generate(std::string_view operator_name, const IO::InputBlock &input,
         const Wavefunction &wf);

//! List available operators
void list_operators();

} // namespace DiracOperator
