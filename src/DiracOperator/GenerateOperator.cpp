#include "GenerateOperator.hpp"
#include "IO/InputBlock.hpp"
#include "fmt/color.hpp"
#include "qip/String.hpp"
#include <iostream>

namespace DiracOperator {

//--------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate(std::string_view operator_name, const IO::InputBlockLegacy &input,
         const Wavefunction &wf) {

  for (const auto &entry : Registry::get().entries()) {
    if (qip::ci_wc_compare(entry.name, operator_name)) {
      return entry.factory(input, wf);
    }
  }

  std::cerr << "\nFAILED to find operator: " << operator_name << "\n";

  const auto name_list = [&]() {
    std::vector<std::string> out;
    for (const auto &entry : Registry::get().entries())
      out.push_back(entry.name);
    return out;
  }();
  std::cout << "Did you mean: "
            << qip::ci_closest_match(operator_name, name_list) << " ?\n ";

  std::cout << "Currently available operators:\n";
  for (const auto &entry : Registry::get().entries())
    std::cout << "  " << entry.name << "\n";
  std::cout << "\n";
  std::cout << "Returning NULL operator (0)\n";
  return std::make_unique<NullOperator>(NullOperator());
}

//--------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate(std::string_view operator_name, const IO::InputBlock &input,
         const Wavefunction &wf) {
  IO::InputBlockLegacy legacy{std::string(operator_name),
                              input.to_ampsci_string()};
  return generate(operator_name, legacy, wf);
}

//--------------------------------------------------------------------
void list_operators() {
  for (const auto &entry : Registry::get().entries())
    fmt::print("  {:10s}  :  {}\n", entry.name, entry.description);
}

} // namespace DiracOperator
