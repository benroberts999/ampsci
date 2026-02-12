#include "GenerateOperator.hpp"

namespace DiracOperator {

//--------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate(std::string_view operator_name, const IO::InputBlock &input,
         const Wavefunction &wf) {

  for (auto &[t_name, t_generator, desc] : operator_list) {
    if (qip::ci_wc_compare(t_name, operator_name)) {
      return t_generator(input, wf);
    }
  }

  std::cerr << "\nFAILED to find operator: " << operator_name
            << " in DiracOperator::generate. Make sure you add it to "
               "operator_list!\n";

  const auto name_list = [&]() {
    std::vector<std::string> out;
    for (const auto &[name, generator, desc] : operator_list) {
      out.push_back(name);
    }
    return out;
  }();
  std::cout << "Did you mean: "
            << *qip::ci_closest_match(operator_name, name_list) << " ?\n ";

  std::cout << "Currently available operators:\n";
  for (const auto &[name, generator, desc] : operator_list) {
    std::cout << "  " << name << "\n";
  }
  std::cout << "\n";
  std::cout << "Returning NULL operator (0)\n";
  return std::make_unique<NullOperator>(NullOperator());
}

//--------------------------------------------------------------------
void list_operators() {
  for (auto &[name, func, description] : operator_list) {
    // std::cout << "  " << name << " " << description << '\n';
    fmt::print("  {:10s}  :  {}\n", name, description);
  }
}

} // namespace DiracOperator