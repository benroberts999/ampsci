#pragma once
#include <type_traits>

namespace qip {

//! Helper template for comparisons. Derive from this to provide !=,>,<=,>=, given == and <.
/*!
If you have a class, C, derive publically from this, with T=C.
Then, just implement '==' and '<' operators, to magically get all of these.

class C : public qip::Comparison<C> {
  //...
  friend bool operator==(const C &lhs, const C &rhs);
  friend bool operator<(const C &lhs, const C &rhs);
};

Can have different types for LHS and RHS, e.g.:

class D : public qip::Comparison<D>, 
                 qip::Comparison<D,C>, 
                qip::Comparison<C,D> {
  //...
  friend bool operator==(const D &lhs, const D &rhs);
  friend bool operator<(const D &lhs, const D &rhs);
  friend bool operator==(const D &lhs, const C &rhs);
  friend bool operator<(const D &lhs, const C &rhs);
  friend bool operator==(const C &lhs, const D &rhs);
  friend bool operator<(const C &lhs, const D &rhs);
};
*/
template <typename T, typename U = T> class Comparison {
  friend bool operator!=(const T &lhs, const U &rhs) { return !(lhs == rhs); }
  friend bool operator>(const T &lhs, const U &rhs) { return rhs < lhs; }
  friend bool operator<=(const T &lhs, const U &rhs) { return !(lhs > rhs); }
  friend bool operator>=(const T &lhs, const U &rhs) { return !(lhs < rhs); }
};

//! Helper template for Arithmetic operations. Derive from this to provide +,-,*,/, given +=, -=, *=, /=
template <typename T> class Arithmetic {
  friend T operator+(T lhs, const T &rhs) { return lhs += rhs; }
  friend T operator-(T lhs, const T &rhs) { return lhs -= rhs; }
  friend T operator*(T lhs, const T &rhs) { return lhs *= rhs; }
  friend T operator/(T lhs, const T &rhs) { return lhs /= rhs; }
};

//! Helper template for Arithmetic operations. Derive from this to provide +,-,*,/, given +=, -=, *=, /=. Works for two different types
template <typename T, typename U> class Arithmetic2 {
  friend T operator+(T lhs, const U &rhs) { return lhs += rhs; }
  friend T operator-(T lhs, const U &rhs) { return lhs -= rhs; }
  friend T operator*(T lhs, const U &rhs) { return lhs *= rhs; }
  friend T operator/(T lhs, const U &rhs) { return lhs /= rhs; }
  friend T operator+(const U &lhs, T rhs) { return rhs += lhs; }
  friend T operator-(const U &lhs, T rhs) { return rhs -= lhs; }
  friend T operator*(const U &lhs, T rhs) { return rhs *= lhs; }
  friend T operator/(const U &lhs, T rhs) { return rhs /= lhs; }
};

} // namespace qip