#pragma once
#include <iostream>
#include <type_traits>

//! Collection of handy tools
namespace qip {

/*!
@brief A light-weight easy-to-use single-file header-only template class for
strong typing.

@details
It works using scoped enums; the enum value (+scope) are used to ensure
uniqueness for a given user-defined strong type. The class template, defined in
the qip namespace, takes an enum class value, and a base type. The base type is
confined to be an arithmetic type (float, double, int, long, etc.) to allow for
operator overloading. The 'using' declaration is optional, but makes things much
easier.

```cpp
using MyType = qip::StrongType<EnumClass::value, BaseType>;
```

The newly-defined Strong Type will, for the most part, behave just like an
instance of the underlying base type would. In particular, the usual arithmetic
operations (+, -, \*, /, +=, \*= etc.) are all defined, and they work with
iostreams. The main difference is that all implicit conversions are banned.

## Examples

Usage is best shown with examples. Consider this dummy problem, where we define
strong types for Energy, Mass, and Velocity.

```cpp
#include "StrongType.hpp"
enum class MechanicsTypes { energy, mass, velocity };
// Use of enum class ensures each StrongType is unique

using Energy   = qip::StrongType<MechanicsTypes::energy, double>;
using Mass     = qip::StrongType<MechanicsTypes::mass, double>;
using Velocity = qip::StrongType<MechanicsTypes::velocity, double>;
```
*/
template <auto enumV, typename BaseT>
struct StrongType {
private:
  static_assert(std::is_arithmetic_v<BaseT>,
                "StrongType only available for arithmetic types");
  static_assert(
    std::is_enum_v<decltype(enumV)>,
    "StrongType must be instantiated with scoped enum (enum class)");
  using StrongT = StrongType<enumV, BaseT>; // type alias

public:
  BaseT v;

  explicit constexpr StrongType(BaseT tv) : v(tv) {}
  explicit constexpr operator BaseT() const { return v; }
  constexpr BaseT &as_base() { return v; }
  [[nodiscard]] constexpr BaseT as_base() const { return v; }

  //! makes 'BaseType' publicly accessible
  using BaseType = BaseT;

  //! Provides operators for regular arithmetic operations
  constexpr StrongT &operator*=(const StrongT &rhs) {
    this->v *= rhs.v;
    return *this;
  }
  friend constexpr StrongT operator*(StrongT lhs, const StrongT &rhs) {
    return lhs *= rhs;
  }
  constexpr StrongT &operator/=(const StrongT &rhs) {
    this->v /= rhs.v;
    return *this;
  }
  friend constexpr StrongT operator/(StrongT lhs, const StrongT &rhs) {
    return lhs /= rhs;
  }
  constexpr StrongT &operator+=(const StrongT &rhs) {
    this->v += rhs.v;
    return *this;
  }
  friend constexpr StrongT operator+(StrongT lhs, const StrongT &rhs) {
    return lhs += rhs;
  }
  constexpr StrongT &operator-=(const StrongT &rhs) {
    this->v -= rhs.v;
    return *this;
  }
  friend constexpr StrongT operator-(StrongT lhs, const StrongT &rhs) {
    return lhs -= rhs;
  }

  //! Provide Base*Strong, Strong*Base oprators - allow scalar multiplication
  constexpr StrongT &operator*=(const BaseT &rhs) {
    this->v *= rhs;
    return *this;
  }
  friend constexpr StrongT operator*(StrongT lhs, const BaseT &rhs) {
    return lhs *= rhs;
  }
  friend constexpr StrongT operator*(const BaseT &lhs, StrongT rhs) {
    return rhs *= lhs;
  }
  //! Provide Strong/Base, but NOT Base/Strong (still scalar multiplication).
  // If StrongT is used for physical units, this will likely not be what you
  // want. In this case, just be explicit. Base/Strong is not scalar
  // multiplication.
  constexpr StrongT &operator/=(const BaseT &rhs) {
    this->v /= rhs;
    return *this;
  }
  friend constexpr StrongT operator/(StrongT lhs, const BaseT &rhs) {
    return lhs /= rhs;
  }

  //! Provides pre/post increment/decrement (++, --) operators
  constexpr StrongT &operator++() {
    ++v;
    return *this;
  }
  constexpr StrongT operator++(int) {
    StrongT result(*this);
    ++(*this);
    return result;
  }
  constexpr StrongT &operator--() {
    --v;
    return *this;
  }
  constexpr StrongT operator--(int) {
    StrongT result(*this);
    --(*this);
    return result;
  }

  //! Provides comparison operators
  friend constexpr bool operator==(const StrongT &lhs, const StrongT &rhs) {
    return lhs.v == rhs.v;
  }
  friend constexpr bool operator!=(const StrongT &lhs, const StrongT &rhs) {
    return !(lhs == rhs);
  }
  friend constexpr bool operator<(const StrongT &lhs, const StrongT &rhs) {
    return lhs.v < rhs.v;
  }
  friend constexpr bool operator>(const StrongT &lhs, const StrongT &rhs) {
    return rhs < lhs;
  }
  friend constexpr bool operator<=(const StrongT &lhs, const StrongT &rhs) {
    return !(rhs < lhs);
  }
  friend constexpr bool operator>=(const StrongT &lhs, const StrongT &rhs) {
    return !(lhs < rhs);
  }

  //! Provides operators for direct comparison w/ BaseT literal (rvalue).
  //! Note: Does not allow comparison with BaseT lvalue
  friend constexpr bool operator==(const StrongT &lhs, const BaseT &&rhs) {
    return lhs.v == rhs;
  }
  friend constexpr bool operator!=(const StrongT &lhs, const BaseT &&rhs) {
    return lhs.v != rhs;
  }
  friend constexpr bool operator<(const StrongT &lhs, const BaseT &&rhs) {
    return lhs.v < rhs;
  }
  friend constexpr bool operator>(const StrongT &lhs, const BaseT &&rhs) {
    return lhs.v > rhs;
  }
  friend constexpr bool operator<=(const StrongT &lhs, const BaseT &&rhs) {
    return lhs.v <= rhs;
  }
  friend constexpr bool operator>=(const StrongT &lhs, const BaseT &&rhs) {
    return lhs.v >= rhs;
  }
  friend constexpr bool operator==(const BaseT &&lhs, const StrongT &rhs) {
    return lhs == rhs.v;
  }
  friend constexpr bool operator!=(const BaseT &&lhs, const StrongT &rhs) {
    return lhs != rhs.v;
  }
  friend constexpr bool operator<(const BaseT &&lhs, const StrongT &rhs) {
    return lhs < rhs.v;
  }
  friend constexpr bool operator>(const BaseT &&lhs, const StrongT &rhs) {
    return lhs > rhs.v;
  }
  friend constexpr bool operator<=(const BaseT &&lhs, const StrongT &rhs) {
    return lhs <= rhs.v;
  }
  friend constexpr bool operator>=(const BaseT &&lhs, const StrongT &rhs) {
    return lhs >= rhs.v;
  }

  //! Provides iostream interface, works as it would for BaseT
  friend std::ostream &operator<<(std::ostream &os, const StrongT &rhs) {
    return os << rhs.v;
  }
  friend std::istream &operator>>(std::istream &is, StrongT &rhs) {
    return is >> rhs.v;
  }
};

} // namespace qip
