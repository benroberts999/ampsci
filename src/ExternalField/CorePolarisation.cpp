#include "CorePolarisation.hpp"
#include "DiagramRPA.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "HF/HartreeFock.hpp"
#include "TDHF.hpp"
#include "TDHFbasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/color.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace ExternalField {

//==============================================================================
std::unique_ptr<CorePolarisation>
make_rpa(const std::string &method, const DiracOperator::TensorOperator *h,
         const HF::HartreeFock *vhf, bool print,
         const std::vector<DiracSpinor> &basis, const std::string &identity) {

  // Parse method for RPA:
  auto rpa_method = ExternalField::ParseMethod(method);

  if (rpa_method == ExternalField::Method::Error) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 148: ");
    fmt::print(
      "RPA method {} not found - check spelling? Defaulting to NO rpa\n",
      method);
    rpa_method = ExternalField::Method::none;
  }
  const auto rpaQ = rpa_method != ExternalField::Method::none;

  // do RPA:
  std::unique_ptr<ExternalField::CorePolarisation> dV{nullptr};
  if (rpaQ && print)
    std::cout << "Including RPA: ";
  if (rpa_method == ExternalField::Method::TDHF) {
    if (print)
      std::cout << "TDHF method\n";
    dV = std::make_unique<ExternalField::TDHF>(h, vhf);
  } else if (rpa_method == ExternalField::Method::basis) {
    if (print)
      std::cout << "TDHF/basis method (" << DiracSpinor::state_config(basis)
                << ")\n";
    dV = std::make_unique<ExternalField::TDHFbasis>(h, vhf, basis);
  } else if (rpa_method == ExternalField::Method::diagram) {
    if (print)
      std::cout << "diagram method (" << DiracSpinor::state_config(basis)
                << ")\n";
    dV = std::make_unique<ExternalField::DiagramRPA>(h, basis, vhf, identity);
  }
  return dV;
}

} // namespace ExternalField
