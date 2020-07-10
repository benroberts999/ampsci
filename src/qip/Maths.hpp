#pragma once
#include <type_traits>

namespace qip {

inline auto comp_abs = [](auto a, auto b) { return std::abs(a) < std::abs(b); };

//******************************************************************************
//! x^n for integer n (n compile-time template parameter)
template <int n, typename T> constexpr auto pow(T x) {
  static_assert(
      std::is_arithmetic_v<T>,
      "In compare(std::vector<T>, std::vector<T>), T must be arithmetic");
  // Returns double for inverse powers, T otherwise
  if constexpr (n < 0) {
    return double(1.0) / pow<-n>(x);
  } else if constexpr (n == 0) {
    (void)x; // 'x' unused in this branch
    return static_cast<T>(1);
    // } else if constexpr (n == 1) {
    //   return x;
  } else {
    return x * pow<n - 1>(x);
  }
}

} // namespace qip
