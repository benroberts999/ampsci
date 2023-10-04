#pragma once
#include "qip/Methods.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <utility>
#include <vector>

namespace qip {

// Uses "no non-const ref" rule. Pass by pointer means modify value

//==============================================================================
//! Merges a number of vectors: {a,b,c},{d,e},{f} -> {a,b,c,d,e,f}
template <typename T, typename... Args>
std::vector<T> merge(std::vector<T> first, const std::vector<T> &second,
                     const Args &...rest) {
  first.insert(first.end(), second.cbegin(), second.cend());
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    return merge(std::move(first), rest...);
  }
}
//==============================================================================
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

//==============================================================================
//! Adds any number of vectors, in place (modifies first vector). Must be of
//! same type. May allocate; will resize first to be size of largest vector.
template <typename T, typename... Args>
void add(std::vector<T> *first, const std::vector<T> &second,
         const Args &...rest) {
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
add(std::vector<T> first, const std::vector<T> &second, const Args &...rest) {
  add(&first, second, rest...);
  return first;
}

// template <typename T, typename... Args>
// [[nodiscard]] std::vector<T> add(std::vector<T> first, const Args &... rest)
// {
//   add(&first, rest...);
//   return first;
// }

//==============================================================================
//! Multiplies any number of (arithmetic) vectors, in place (modifies first
//! vector). Must be of same type. May allocate; will resize first to be size of
//! largest vector.
template <typename T, typename... Args>
void multiply(std::vector<T> *first, const std::vector<T> &second,
              const Args &...rest) {

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
                                      const Args &...rest) {
  multiply(&first, second, rest...);
  return first;
}

//==============================================================================
//! Composes any number of vectors, in place (modifies first vector), using the
//! provided function. Must be of same type. May allocate; will resize first to
//! be size of largest vector.
//! e.g., qip::compose(std::plus{}, &vo, v2, v3); same as qip::add(&vo, v2, v3)
template <typename F, typename T, typename... Args>
void compose(const F &func, std::vector<T> *first, const std::vector<T> &second,
             const Args &...rest) {
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
                                     const Args &...rest) {
  compose(func, &first, second, rest...);
  return first;
}

//==============================================================================
//! In-place scalar multiplication of std::vector - types must match
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

//==============================================================================
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

//! Produces a Log-Linear distributed range of values between
//! [first,last] with number steps. number must be at least 3. first+last are
//! guarenteed to be the first and last points in the range. T must be floating
//! point. Range is roughly logarithmic for values below 'b', and linear for
//! values above b. Not tested for negative values.
template <typename T, typename N>
std::vector<T> loglinear_range(T first, T last, T b, N number) {
  static_assert(std::is_floating_point_v<T>,
                "In loglinear_range(T, T, T, N), T must be floating point");
  static_assert(std::is_integral_v<N>,
                "In loglinear_range(T, T, T, N), N must be integral");
  assert(number >= 3);
  assert(first > 0.0 && last > 0.0 && b > 0.0);

  std::vector<T> range;
  range.reserve(static_cast<std::size_t>(number));

  const auto du = (last - first + b * std::log(last / first)) /
                  (static_cast<T>(number - 1));

  auto next_r = [b](auto u1, auto r_guess) {
    // Solve eq. u = r + b ln(r) to find r
    // Use Newtons method
    // => f(r) = r + b ln(r) - u = 0
    // dfdr = b(1/r + 1/(b+r))
    const auto f_u = [b, u1](double tr) { return tr + b * std::log(tr) - u1; };
    const auto dr = 0.1 * r_guess;
    const auto delta_targ = r_guess * 1.0e-18;
    const auto [ri, delta_r] = qip::Newtons(f_u, r_guess, delta_targ, dr, 30);
    return ri;
  };

  auto u = first + b * std::log(first);
  range.push_back(first);
  for (N i = 1; i < number - 1; ++i) {
    u += du;
    range.push_back(next_r(u, range.back()));
  }
  range.push_back(last);

  return range;
}

//==============================================================================
//! first[i]*...rest[i]  --  used to allow inner_product
template <typename T, typename... Args>
constexpr auto multiply_at(std::size_t i, const T &first, const Args &...rest) {
  // assert(first.size() < i);
  if constexpr (sizeof...(rest) == 0) {
    return first[i];
  } else {
    return first[i] * multiply_at(i, rest...);
  }
}

//! Variadic inner product (v1,v2,...vn) : sum_i v1[i]*v2[i]*...vn[i]
template <typename T, typename... Args>
constexpr auto inner_product(const T &first, const Args &...rest) {
  auto res = multiply_at(0, first, rest...);
  for (std::size_t i = 1; i < first.size(); ++i) {
    res += multiply_at(i, first, rest...);
  }
  return res;
}

template <typename T, typename... Args>
auto inner_product_sub(std::size_t p0, std::size_t pinf, const T &first,
                       const Args &...rest) {
  auto res = decltype(first[0])(0);
  for (std::size_t i = p0; i < pinf; ++i) {
    res += multiply_at(i, first, rest...);
  }
  return res;
}

//==============================================================================
template <typename F, typename T> T apply_to(const F &func, T list) {
  for (auto &l : list) {
    l = func(l);
  }
  return list;
}

//==============================================================================
//! Creates a subset of a vector, whose element match condition. By copy.
//! condition must have function signature: bool condition(T)
template <typename T, typename Func>
std::vector<T> select_if(const std::vector<T> &in, Func condition) {
  std::vector<T> out;
  for (const auto &el : in) {
    if (condition(el)) {
      out.push_back(el);
    }
  }
  return out;
}

//! Inserts elements from in into inout, if condition is met
template <typename T, typename Func>
void insert_into_if(const std::vector<T> &in, std::vector<T> *inout,
                    Func condition) {
  for (const auto &el : in) {
    if (condition(el)) {
      inout->push_back(el);
    }
  }
}

//==============================================================================
//==============================================================================
//! namespace qip::overloads provides operator overloads for std::vector
namespace overloads {

//! Provide addition of two vectors:
template <typename T>
std::vector<T> &operator+=(std::vector<T> &a, const std::vector<T> &b) {
  // The following allows this to work for types that can't be default constructed
  const auto smaller_size = std::min(a.size(), b.size());
  if (a.size() < b.size())
    a.insert(a.end(), b.begin() + long(a.size()), b.end());
  for (auto i = 0ul; i < smaller_size; ++i) {
    a[i] += b[i];
  }
  return a;
}
template <typename T>
std::vector<T> operator+(std::vector<T> a, const std::vector<T> &b) {
  return a += b;
}
// and subtraction
template <typename T>
std::vector<T> &operator-=(std::vector<T> &a, const std::vector<T> &b) {
  const auto size = std::max(a.size(), b.size());
  a.resize(size);
  for (auto i = 0ul; i < b.size(); ++i) {
    a[i] -= b[i];
  }
  return a;
}
template <typename T>
std::vector<T> operator-(std::vector<T> a, const std::vector<T> &b) {
  return a -= b;
}

// Provide scalar multiplication
template <typename T, typename U>
std::vector<T> &operator*=(std::vector<T> &v, U x) {
  if (x != U{1}) {
    for (auto &v_i : v) {
      v_i *= x;
    }
  }
  return v;
}
template <typename T, typename U>
std::vector<T> operator*(std::vector<T> v, U x) {
  return v *= x;
}
template <typename T, typename U>
std::vector<T> operator*(U x, std::vector<T> v) {
  return v *= x;
}

// Provide scalar devision
template <typename T, typename U>
std::vector<T> &operator/=(std::vector<T> &v, U x) {
  if (x != U{1}) {
    for (auto &v_i : v) {
      v_i /= x;
    }
  }
  return v;
}
template <typename T, typename U>
std::vector<T> operator/(std::vector<T> v, U x) {
  return v /= x;
}

//! In-place element-wise multiplication: a*=b => a_i := a_i * b_i
template <typename T>
std::vector<T> &operator*=(std::vector<T> &a, const std::vector<T> &b) {
  const auto min_size = std::min(a.size(), b.size());
  for (std::size_t i = 0; i < min_size; ++i) {
    a[i] *= b[i];
  }
  for (std::size_t i = min_size; i < a.size(); ++i) {
    a[i] = 0.0;
  }
  if (a.size() < b.size())
    a.resize(b.size());
  return a;
}
//! Element-wise multiplication: c=a*b => c_i = a_i * b_i
template <typename T>
std::vector<T> operator*(std::vector<T> a, const std::vector<T> &b) {
  return a *= b;
}

//! In-place element-wise division: a/=b => a_i := a_i / b_i
template <typename T>
std::vector<T> &operator/=(std::vector<T> &a, const std::vector<T> &b) {
  assert(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    a[i] /= b[i];
  }
  return a;
}
//! Element-wise multiplication: c=a/b => c_i = a_i / b_i
template <typename T>
std::vector<T> operator/(std::vector<T> a, const std::vector<T> &b) {
  return a /= b;
}

} // namespace overloads

} // namespace qip
