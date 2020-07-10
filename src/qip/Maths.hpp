#pragma once
#include <cmath>
#include <type_traits>

namespace qip {

inline auto comp_abs = [](auto a, auto b) { return std::abs(a) < std::abs(b); };

//******************************************************************************
template <typename T, typename... Args> T max_abs(T first, Args... rest) {
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    const auto max_rest = max_abs(rest...);
    if (std::abs(first) >= std::abs(max_rest))
      return first;
    return max_rest;
  }
}
template <typename T, typename... Args> T min_abs(T first, Args... rest) {
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    const auto min_rest = min_abs(rest...);
    if (std::abs(first) <= std::abs(min_rest))
      return first;
    return min_rest;
  }
}

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
