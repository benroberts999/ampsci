#pragma once
#include "IO/InputBlock.hpp"
#include "Operators/Ek.hpp"
#include "Operators/hfs.hpp"
#include "Operators/M1.hpp"
#include "Operators/PNC.hpp"
#include "Operators/QED.hpp"
#include "Operators/RadialF.hpp"
#include "Operators/p.hpp"
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
    operator_list{
        {"E1", &generate_E1}, {"E1v", &generate_E1v}, {"E2", &generate_E2},
        {"Ek", &generate_Ek}, {"M1", &generate_M1},   {"hfs", &generate_hfs},
        {"r", &generate_r},   {"pnc", &generate_pnc}, {"Vrad", &generate_Vrad},
        {"dr", &generate_dr}, {"p", &generate_p}};

//------------------------------------------------------------------------------

//! Returns a unique_ptr (polymorphic) to the requested operator, with given
//! properties
inline std::unique_ptr<DiracOperator::TensorOperator>
generate(std::string_view operator_name, const IO::InputBlock &input,
         const Wavefunction &wf) {

  for (auto &[t_name, t_generator] : operator_list) {
    if (qip::ci_wc_compare(t_name, operator_name)) {
      return t_generator(input, wf);
    }
  }

  std::cerr << "\nFAILED to find operator: " << operator_name
            << " in DiracOperator::generate. Make sure you add it to "
               "operator_list!\n";

  const auto name_list = [&]() {
    std::vector<std::string> out;
    for (const auto &[name, generator] : operator_list) {
      out.push_back(name);
    }
    return out;
  }();
  std::cout << "Did you mean: "
            << *qip::ci_closest_match(operator_name, name_list) << " ?\n ";

  std::cout << "Currently available operators:\n";
  for (const auto &[name, generator] : operator_list) {
    std::cout << "  " << name << "\n";
  }
  std::cout << "\n";
  std::cout << "Returning NULL operator (0)\n";
  return std::make_unique<NullOperator>(NullOperator());
}

inline void list_operators() {
  for (auto &[name, func] : operator_list) {
    std::cout << "  " << name << '\n';
  }
}

} // namespace DiracOperator
