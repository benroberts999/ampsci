#pragma once
#include <cassert>
#include <iomanip>
#include <iostream> //stringstream?
#include <string_view>

namespace qip {

template <typename T>
bool check_value(std::ostream *ostr, const std::string_view &name,
                 const T value, const T expected, const T tollerance) {
  const bool passedQ = std::abs(value - expected) <= tollerance;
  if (passedQ) {
    *ostr << "passed ";
  } else {
    *ostr << "FAILED ";
  }

  auto val_pr = (expected == 0.0) ? std::setprecision(1) : std::setprecision(5);

  *ostr << std::setw(36) << std::left << name << "| " << std::right
        << std::setw(8) << std::scientific << val_pr << value;
  if ((expected != 0.0))
    *ostr << "/" << expected;
  *ostr << " [" << std::setprecision(0) << tollerance << "]\n";

  return passedQ;
}

} // namespace qip
