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

template <typename T, typename... Args> T max(T first, Args... rest) {
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    const auto max_rest = max(rest...);
    if (first >= max_rest)
      return first;
    return max_rest;
  }
}
template <typename T, typename... Args> T min(T first, Args... rest) {
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    const auto min_rest = min(rest...);
    if (first <= min_rest)
      return first;
    return min_rest;
  }
}

template <typename T, typename... Args>
T max_difference(T first, Args... rest) {
  return max(first, rest...) - min(first, rest...);
}

//******************************************************************************
//! x^n for integer n (n compile-time template parameter)
template <int n, typename T> constexpr auto pow(T x) {
  static_assert(std::is_arithmetic_v<T>, "In pow(T x), T must be arithmetic");
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

//! x^n for integer n (runtime n)
template <typename T> constexpr T pow(T x, int n) {
  static_assert(std::is_floating_point_v<T>,
                "In pow(T x, int n), T must be floating point");
  if (n < 0)
    return static_cast<T>(1) / pow<T>(x, -n);
  T result = static_cast<T>(1);
  for (int i = 0; i < n; ++i) {
    result *= x;
  }
  return result;
}

//******************************************************************************
//! Returns sign of value
template <typename T> constexpr int sign(T value) {
  static_assert(std::is_arithmetic_v<T>,
                "In sign(T value), T must be arithmetic");
  return (T(0) < value) - (value < T(0));
}

//! Clips value to between -max <= value <= max
//! clip(x,max) : |x| > max, ret max; if |x|<-max, -max; else x
template <typename T> constexpr T clip(T value, T max_abs) {
  static_assert(std::is_arithmetic_v<T>,
                "In clip(T value, T max_abs), T must be arithmetic");
  if (value > max_abs)
    return max_abs;
  if (value < -max_abs)
    return -max_abs;
  return value;
}

//! Sets values |v|<min to zero; if |v|>=min, returns v
template <typename T> constexpr T chop(T value, T min_abs) {
  static_assert(std::is_arithmetic_v<T>,
                "In clip(T value, T max_abs), T must be arithmetic");
  if (std::abs(value) < min_abs)
    return static_cast<T>(0);
  return value;
}

} // namespace qip
