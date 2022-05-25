#pragma once
#include <cassert>
#include <iomanip>
#include <iostream> //stringstream?
#include <string_view>

namespace qip {

//==============================================================================
//! Checks if value==expected, within +/-tollerance. Writes report to ostr.
template <typename T>
bool check_value(std::ostream *ostr, const std::string_view &name,
                 const T value, const T expected, const T tollerance) {
  static_assert(
      std::is_arithmetic_v<T>,
      "In check_value(os *o, sv s, T x, T y, T z), T must be arithmetic");
  const bool passedQ = std::abs(value - expected) <= tollerance;
  if (passedQ) {
    *ostr << "passed ";
  } else {
    *ostr << "FAILED ";
  }

  const auto val_pr =
      (expected == T{0}) ? std::setprecision(1) : std::setprecision(5);

  *ostr << std::setw(36) << std::left << name << "| " << std::right
        << std::setw(8) << std::scientific << val_pr << value;
  if ((expected != T{0}))
    *ostr << "/" << expected;
  *ostr << " [" << std::setprecision(0) << tollerance << "]\n";

  return passedQ;
}

//! overload which prints to cout
template <typename T>
bool check(const std::string_view &name, const T value, const T expected) {
  return check(&std::cout, name, value, expected);
}

template <typename T>
bool check(std::ostream *ostr, const std::string_view &name, const T value,
           const T expected) {
  const bool passedQ = value == expected;
  if (passedQ) {
    *ostr << "passed ";
  } else {
    *ostr << "FAILED ";
  }

  *ostr << std::setw(36) << std::left << name << "|\n";

  return passedQ;
}

inline bool check(const std::string_view &name, const bool condition) {
  return check<bool>(&std::cout, name, condition, true);
}

inline bool check(std::ostream *ostr, const std::string_view &name,
                  const bool condition) {
  return check<bool>(ostr, name, condition, true);
}

} // namespace qip
