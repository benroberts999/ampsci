#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <utility>
#include <vector>

namespace qip {

// Uses "no non-const ref" rule. Pass by pointer means modify value

//******************************************************************************
//! Directly compare two arithmetic vectors of the same type and length.
//! Returns pair {delta, itr} where delta = |max|{first - second}, itr is
//! iterator to position in first vector where the maximum delta occured.
//! Note: Maximum is by magnitude, but delta is signed as (first-second)
template <typename T>
auto compare(const std::vector<T> &first, const std::vector<T> &second) {
  static_assert(
      std::is_arithmetic_v<T>,
      "In compare(std::vector<T>, std::vector<T>), T must be arithmetic");
  assert(first.size() == second.size()); //?

  auto it1 = first.cbegin();
  auto it2 = second.cbegin();
  auto largest_delta = static_cast<T>(0);
  auto largest_at = first.cend();
  for (; it1 != first.cend(); ++it1, ++it2) {
    const auto delta = *it1 - *it2;
    if (std::abs(delta) > std::abs(largest_delta)) {
      largest_delta = delta;
      largest_at = it1;
    }
  }
  return std::make_pair(largest_delta, largest_at);
}

//! Compares two vectors of the same length, according to the rule given by
//! func. Returns pair {delta, itr} where delta = max{|func(first,second)|}, itr
//! is iterator to position in first vector where the maximum delta occured.
//! Note: Maximum is by magnitude.
template <typename T, typename U, typename Func>
auto compare(const std::vector<T> &first, const std::vector<U> &second,
             Func &func) {
  assert(first.size() == second.size()); //?

  auto it1 = first.cbegin();
  auto it2 = second.cbegin();
  auto largest_eps = 0.0; // double always?
  auto largest_at = first.cend();
  for (; it1 != first.cend(); ++it1, ++it2) {
    const auto eps = func(*it1, *it2);
    if (std::abs(eps) > std::abs(largest_eps)) {
      largest_eps = eps;
      largest_at = it1;
    }
  }
  return std::make_pair(largest_eps, largest_at);
}

//! Compares values of two arithmetic vectors of the same type and
//! length, relative to second value.
//! Returns pair {eps, itr} where eps = |max|{(first - second)/second},
//! itr is iterator to position in first vector where the maximum eps occured.
//! Note: Maximum is by magnitude, but eps is signed as (first-second)/second
template <typename T>
auto compare_eps(const std::vector<T> &first, const std::vector<T> &second) {
  static_assert(
      std::is_floating_point_v<T>,
      "In compare_eps(std::vector<T>, std::vector<T>), T must be floating pt");
  assert(first.size() == second.size()); //?

  auto it1 = first.cbegin();
  auto it2 = second.cbegin();
  auto largest_eps = static_cast<T>(0);
  auto largest_at = first.cend();
  for (; it1 != first.cend(); ++it1, ++it2) {
    const auto eps = (*it1 - *it2) / (*it2); // may be +/-inf, fine
    if (std::abs(eps) > std::abs(largest_eps)) {
      largest_eps = eps;
      largest_at = it1;
    }
  }
  return std::make_pair(largest_eps, largest_at);
}

//******************************************************************************
//! Adds any number of vectors, in place (modifies first vector). Must be of
//! same type. May allocate; will resize first to be size of largest vector.
template <typename T, typename... Args>
void add(std::vector<T> *first, const std::vector<T> &second,
         const Args &... rest) {
  // re-sizes vector:
  if (first->size() < second.size())
    first->resize(second.size());

  for (std::size_t i = 0; i < second.size(); ++i) {
    (*first)[i] += second[i];
  }
  if constexpr (sizeof...(rest) != 0)
    add(first, rest...);
}

//! Adds any number of vectors. Must be of same type. Size of returned vector
//! will be size of largest vector; may allocate more than once if vectors are
//! not all of same size.
template <typename T, typename... Args>
[[nodiscard]] std::vector<T>
add(std::vector<T> first, const std::vector<T> &second, const Args &... rest) {
  add(&first, second, rest...);
  return first;
}

// template <typename T, typename... Args>
// [[nodiscard]] std::vector<T> add(std::vector<T> first, const Args &... rest)
// {
//   add(&first, rest...);
//   return first;
// }

//******************************************************************************
//! Multiplies any number of (arithmetic) vectors, in place (modifies first
//! vector). Must be of same type. May allocate; will resize first to be size of
//! largest vector.
template <typename T, typename... Args>
void multiply(std::vector<T> *first, const std::vector<T> &second,
              const Args &... rest) {

  static_assert(std::is_arithmetic_v<T>,
                "In multiply(std::vector<T>...) : T must be arithmetic");

  // {1,1,1,1}*{1,1} -> {1,1,0,0}; re-sizes vector
  const auto size = std::min(first->size(), second.size());
  if (first->size() < second.size())
    first->resize(second.size());

  for (std::size_t i = 0; i < size; ++i) {
    (*first)[i] *= second[i];
  }
  for (std::size_t i = size; i < first->size(); ++i) {
    (*first)[i] = static_cast<T>(0);
  }
  if constexpr (sizeof...(rest) != 0)
    multiply(first, rest...);
}

//! Multiplies any number of vectors. Must be of same type. Size of returned
//! vector will be size of largest vector; may allocate more than once if
//! vectors are not all of same size.
template <typename T, typename... Args>
[[nodiscard]] std::vector<T> multiply(std::vector<T> first,
                                      const std::vector<T> &second,
                                      const Args &... rest) {
  multiply(&first, second, rest...);
  return first;
}

//******************************************************************************
//! Composes any number of vectors, in place (modifies first vector), using the
//! provided function. Must be of same type. May allocate; will resize first to
//! be size of largest vector.
//! e.g., qip::compose(std::plus{}, &vo, v2, v3); same as qip::add(&vo, v2, v3)
template <typename F, typename T, typename... Args>
void compose(const F &func, std::vector<T> *first, const std::vector<T> &second,
             const Args &... rest) {
  // XXX Comiple=-time constraints on f! Must be T(T,T)!

  if (first->size() < second.size())
    first->resize(second.size());

  for (std::size_t i = 0; i < second.size(); ++i) {
    (*first)[i] = func((*first)[i], second[i]);
  }
  for (std::size_t i = second.size(); i < first->size(); ++i) {
    (*first)[i] = func((*first)[i], static_cast<T>(0));
  }
  if constexpr (sizeof...(rest) != 0)
    compose(func, first, rest...);
}

//! Composes any number of vectors. Must be of same type. Size of returned
//! vector will be size of largest vector; may allocate more than once if
//! vectors are not all of same size.
template <typename F, typename T, typename... Args>
[[nodiscard]] std::vector<T> compose(const F &func, std::vector<T> first,
                                     const std::vector<T> &second,
                                     const Args &... rest) {
  compose(func, &first, second, rest...);
  return first;
}

//******************************************************************************
//! In-pace scalar multiplication of std::vector - types must match
template <typename T> void scale(std::vector<T> *vec, T x) {
  static_assert(std::is_arithmetic_v<T>,
                "In scale(std::vector<T>, T) : T must be arithmetic");
  for (auto &v : *vec)
    v *= x;
}

//! Scalar multiplication of std::vector - types must match
template <typename T>
[[nodiscard]] std::vector<T> scale(std::vector<T> vec, T x) {
  scale(&vec, x);
  return vec;
}

//******************************************************************************
//! Produces a uniformly*(see below) distributed range of values between
//! [first,last] with number steps. number must be at least 2. first+last are
//! guarenteed to be the first and last points in the range. For integral T,
//! range will not be perfectly uniform, due to [first, ..., last] guarentee and
//! rounding; also in this case same value may appear more than once if too-many
//! steps are requested.
template <typename T, typename N>
std::vector<T> uniform_range(T first, T last, N number) {
  static_assert(std::is_arithmetic_v<T>,
                "In uniform_range(T, T, N), T must be arithmetic");
  static_assert(std::is_integral_v<N>,
                "In uniform_range(T, T, N), N must be integral");
  assert(number >= 2);

  std::vector<T> range;
  range.reserve(static_cast<std::size_t>(number));
  const auto interval = static_cast<double>(last - first);
  range.push_back(first); // guarentee first is first
  for (N i = 1; i < number - 1; ++i) {
    const auto eps = static_cast<double>(i) / static_cast<double>(number - 1);
    const auto value = static_cast<double>(first) + (eps * interval);
    range.push_back(static_cast<T>(value));
  }
  range.push_back(last); // guarentee last is last
  return range;
}

//! Produces a logarithmicly*(see below) distributed range of values between
//! [first,last] with number steps. number must be at least 2. first+last are
//! guarenteed to be the first and last points in the range. For integral T,
//! range will not be perfectly logarithmic, due to [first, ..., last] guarentee
//! and rounding; also in this case same value may appear more than once if
//! too-many steps are requested.
template <typename T, typename N>
std::vector<T> logarithmic_range(T first, T last, N number) {
  static_assert(std::is_arithmetic_v<T>,
                "In logarithmic_range(T, T, N), T must be arithmetic");
  static_assert(std::is_integral_v<N>,
                "In logarithmic_range(T, T, N), N must be integral");
  assert(number >= 2);

  std::vector<T> range;
  range.reserve(static_cast<std::size_t>(number));
  const auto log_ratio =
      std::log(static_cast<double>(last) / static_cast<double>(first));
  range.push_back(first);
  for (N i = 1; i < number - 1; ++i) {
    const auto eps = static_cast<double>(i) / static_cast<double>(number - 1);
    const auto value = static_cast<double>(first) * std::exp(log_ratio * eps);
    range.push_back(static_cast<T>(value));
  }
  range.push_back(last);
  return range;
}

} // namespace qip
