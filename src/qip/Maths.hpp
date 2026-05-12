#pragma once
#include "Vector.hpp"
#include <cmath>
#include <type_traits>

namespace qip {

/*!
  @brief Function object for comparisons by absolute value.

  @details
  Compares |lhs| < |rhs| using std::abs. Works like std::less.
  T must be arithmetic.
*/
struct less_abs {
  template <class T>
  constexpr bool operator()(const T &lhs, const T &rhs) const {
    static_assert(std::is_arithmetic_v<T>,
                  "In less_abs<T>(), T must be arithmetic");
    return std::abs(lhs) < std::abs(rhs);
  }
};

//==============================================================================

//! Returns the maximum of any number of parameters (variadic).
template <typename T, typename... Args>
T max(T first, Args... rest) {
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    const auto max_rest = max(rest...);
    if (first >= max_rest)
      return first;
    return max_rest;
  }
}

//! Returns the minimum of any number of parameters (variadic).
template <typename T, typename... Args>
T min(T first, Args... rest) {
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    const auto min_rest = min(rest...);
    if (first <= min_rest)
      return first;
    return min_rest;
  }
}

//! Returns the value with maximum absolute value of any number of parameters (variadic).
template <typename T, typename... Args>
T max_abs(T first, Args... rest) {
  static_assert(std::is_arithmetic_v<T>,
                "In max_abs<T>(), T must be arithmetic");
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    const auto max_rest = max_abs(rest...);
    if (std::abs(first) >= std::abs(max_rest))
      return first;
    return max_rest;
  }
}

//! Returns the value with minimum absolute value of any number of parameters (variadic).
template <typename T, typename... Args>
T min_abs(T first, Args... rest) {
  static_assert(std::is_arithmetic_v<T>,
                "In min_abs<T>(), T must be arithmetic");
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    const auto min_rest = min_abs(rest...);
    if (std::abs(first) <= std::abs(min_rest))
      return first;
    return min_rest;
  }
}

//! Returns max - min over any number of arguments (variadic).
template <typename T, typename... Args>
T max_difference(T first, Args... rest) {
  static_assert(std::is_arithmetic_v<T>,
                "In max_difference<T>(), T must be arithmetic");
  return max(first, rest...) - min(first, rest...);
}

//==============================================================================
/*!
  @brief x^n for compile-time integer n, x any arithmetic type.

  @details
  Returns double for negative n, T otherwise.
*/
template <int n, typename T>
constexpr auto pow(T x) {
  using namespace qip::overloads;
  // Returns double for inverse powers, T otherwise
  if constexpr (n < 0) {
    return double(1.0) / pow<-n>(x);
  } else if constexpr (n == 0) {
    (void)x; // 'x' unused in this branch
    return static_cast<T>(1);
  } else if constexpr (n == 1) {
    return x;
  } else {
    return x * pow<n - 1>(x);
  }
}

//! x^n for runtime integer n; T must be floating point.
template <typename T>
constexpr T pow(T x, int n) {
  static_assert(std::is_floating_point_v<T>,
                "In pow(T x, int n), T must be foating point");
  if (n < 0) {
    return static_cast<T>(1) / pow<T>(x, -n);
  }
  // T result = static_cast<T>(1);
  T result{1};
  for (int i = 0; i < n; ++i) {
    result *= x;
  }
  return result;
}

//==============================================================================

//! Factorial x! - takes integer, returns double.
template <typename T>
constexpr double factorial(T x) {
  static_assert(std::is_integral_v<T>,
                "In factorial(T), T must be an integral type");
  return (x <= 1) ? 1.0 : double(x) * factorial<T>(x - 1);
}

//! Double factorial x!! - takes integer, returns double.
template <typename T>
constexpr double double_factorial(T x) {
  static_assert(std::is_integral_v<T>,
                "double_factorial(T): T must be an integral type");

  return (x <= 1) ? 1.0 : double(x) * double_factorial<T>(x - 2);
}

//==============================================================================

//! Returns the sign of value. Note: sign(0) == 0.
template <typename T>
constexpr int sign(T value) {
  static_assert(std::is_arithmetic_v<T>,
                "In sign(T value), T must be arithmetic");
  return (T(0) < value) - (value < T(0));
}

//! Clamps value to [-max_abs, max_abs] - i.e., using abs.
template <typename T>
constexpr T clamp_abs(T value, T max_abs) {
  static_assert(std::is_arithmetic_v<T>,
                "In clamp_abs(T value, T max_abs), T must be arithmetic");
  if (value > max_abs)
    return max_abs;
  if (value < -max_abs)
    return -max_abs;
  return value;
}

//! Sets value to zero if |value| < min_abs, otherwise returns value unchanged.
template <typename T>
constexpr T chop(T value, T min_abs) {
  static_assert(std::is_arithmetic_v<T>,
                "In cjop(T value, T min_abs), T must be arithmetic");
  if (std::abs(value) < min_abs)
    return static_cast<T>(0);
  return value;
}

} // namespace qip
