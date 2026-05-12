#pragma once
#include <type_traits>

namespace qip {

/*!
  @brief Helper template that provides !=, >, <=, >= given == and <.

  @details
  Derive publicly from this (with T = your class) and implement == and <;
  the remaining comparison operators are provided automatically.

  @code{.cpp}
  class C : public qip::Comparison<C> {
    friend bool operator==(const C &lhs, const C &rhs);
    friend bool operator<(const C &lhs, const C &rhs);
  };
  @endcode

  For mixed-type comparisons, inherit multiple specialisations:

  @code{.cpp}
  class D : public qip::Comparison<D>,
                   qip::Comparison<D,C>,
                   qip::Comparison<C,D> {
    friend bool operator==(const D &lhs, const D &rhs);
    friend bool operator<(const D &lhs, const D &rhs);
    friend bool operator==(const D &lhs, const C &rhs);
    friend bool operator<(const D &lhs, const C &rhs);
    friend bool operator==(const C &lhs, const D &rhs);
    friend bool operator<(const C &lhs, const D &rhs);
  };
  @endcode
*/
template <typename T, typename U = T>
class Comparison {
  friend bool operator!=(const T &lhs, const U &rhs) { return !(lhs == rhs); }
  friend bool operator>(const T &lhs, const U &rhs) { return rhs < lhs; }
  friend bool operator<=(const T &lhs, const U &rhs) { return !(lhs > rhs); }
  friend bool operator>=(const T &lhs, const U &rhs) { return !(lhs < rhs); }
};

/*!
  @brief Helper template that provides +, -, *, / given +=, -=, *=, /=.

  @details
  Derive publicly from this (with T = your class) and implement the
  compound-assignment operators; the binary operators are provided automatically.
*/
template <typename T>
class Arithmetic {
  friend T operator+(T lhs, const T &rhs) { return lhs += rhs; }
  friend T operator-(T lhs, const T &rhs) { return lhs -= rhs; }
  friend T operator*(T lhs, const T &rhs) { return lhs *= rhs; }
  friend T operator/(T lhs, const T &rhs) { return lhs /= rhs; }
};

/*!
  @brief Like @ref Arithmetic, but for two different types (T op U).

  @details
  Provides T+U, T-U, T*U, T/U (and the symmetric U+T, U*T) given the
  compound-assignment operators on T.
*/
template <typename T, typename U>
class Arithmetic2 {
  friend T operator+(T lhs, const U &rhs) { return lhs += rhs; }
  friend T operator-(T lhs, const U &rhs) { return lhs -= rhs; }
  friend T operator*(T lhs, const U &rhs) { return lhs *= rhs; }
  friend T operator/(T lhs, const U &rhs) { return lhs /= rhs; }
  // Also define these symmetric ones
  friend T operator+(const U &lhs, T rhs) { return rhs += lhs; }
  friend T operator*(const U &lhs, T rhs) { return rhs *= lhs; }
};

} // namespace qip
