#include "qip/Template.hpp"
#include "catch2/catch.hpp"
#include <cassert>
#include <iostream>

//-----------------------
class A : public qip::Comparison<A>, qip::Arithmetic<A> {
public:
  A(int x, int y) : a(x), b(y){};
  int a{};
  int b{};
  friend bool operator==(const A &lhs, const A &rhs) { return lhs.a == rhs.a; }
  friend bool operator<(const A &lhs, const A &rhs) { return lhs.a < rhs.a; }

  A &operator+=(const A &rhs) {
    this->a += rhs.a;
    return *this;
  }
  A &operator-=(const A &rhs) {
    this->a -= rhs.a;
    return *this;
  }
  A &operator*=(const A &rhs) {
    this->a *= rhs.a;
    return *this;
  }
};

//-----------------------
class B : qip::Comparison<B>, qip::Comparison<A, B>, qip::Comparison<B, A> {
public:
  B(int x) : a(x){};
  int a{};
  friend bool operator==(const B &lhs, const B &rhs) { return lhs.a == rhs.a; }
  friend bool operator<(const B &lhs, const B &rhs) { return lhs.a < rhs.a; }

  friend bool operator==(const B &lhs, const A &rhs) { return lhs.a == rhs.a; }
  friend bool operator<(const B &lhs, const A &rhs) { return lhs.a < rhs.a; }

  friend bool operator==(const A &lhs, const B &rhs) { return lhs.a == rhs.a; }
  friend bool operator<(const A &lhs, const B &rhs) { return lhs.a < rhs.a; }
};

//-----------------------
class C : qip::Comparison<C>, qip::Arithmetic<C>, qip::Arithmetic2<C, int> {
public:
  explicit C(int x) : a(x){};
  int a{};

  friend bool operator==(const C &lhs, const C &rhs) { return lhs.a == rhs.a; }

  C &operator+=(const C &rhs) {
    this->a += rhs.a;
    return *this;
  }
  C &operator-=(const C &rhs) {
    this->a -= rhs.a;
    return *this;
  }
  C &operator+=(int rhs) {
    this->a += rhs;
    return *this;
  }
  C &operator-=(int rhs) {
    this->a -= rhs;
    return *this;
  }

  C &operator*=(int rhs) {
    this->a *= rhs;
    return *this;
  }
  C &operator/=(int rhs) {
    this->a /= rhs;
    return *this;
  }
};

//==============================================================================
TEST_CASE("qip::Template", "[qip][Template][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "qip::Template\n";

  A a1{6, 1};
  A a2{6, 2};
  A a3{12, 3};

  REQUIRE(a1 == a2);
  REQUIRE_FALSE(a1 != a2);
  REQUIRE(a1 >= a2);
  REQUIRE(a1 <= a2);
  REQUIRE_FALSE(a1 > a2);
  REQUIRE_FALSE(a1 < a2);
  REQUIRE(a1 < a3);
  REQUIRE(a1 <= a3);
  REQUIRE_FALSE(a1 > a3);
  REQUIRE_FALSE(a1 >= a3);

  B b1{6};
  B b2{6};
  B b3{12};

  REQUIRE(b1 == b2);
  REQUIRE_FALSE(b1 != b2);
  REQUIRE(b1 >= b2);
  REQUIRE(b1 <= b2);
  REQUIRE_FALSE(b1 > b2);
  REQUIRE_FALSE(b1 < b2);
  REQUIRE(b1 < b3);
  REQUIRE(b1 <= b3);
  REQUIRE_FALSE(b1 > b3);
  REQUIRE_FALSE(b1 >= b3);

  REQUIRE(b1 == a2);
  REQUIRE_FALSE(b1 != a2);
  REQUIRE(b1 >= a2);
  REQUIRE(b1 <= a2);
  REQUIRE_FALSE(b1 > a2);
  REQUIRE_FALSE(b1 < a2);
  REQUIRE(b1 < a3);
  REQUIRE(b1 <= a3);
  REQUIRE_FALSE(b1 > a3);
  REQUIRE_FALSE(b1 >= a3);

  REQUIRE(a1 == b2);
  REQUIRE_FALSE(a1 != b2);
  REQUIRE(a1 >= b2);
  REQUIRE(a1 <= b2);
  REQUIRE_FALSE(a1 > b2);
  REQUIRE_FALSE(a1 < b2);
  REQUIRE(a1 < b3);
  REQUIRE(a1 <= b3);
  REQUIRE_FALSE(a1 > b3);
  REQUIRE_FALSE(a1 >= b3);

  REQUIRE(a1 + a2 == a3);
  REQUIRE(a1 - a2 == A{0, 1});
  REQUIRE(a1 * a2 == A{36, 1});
  // Purposly didn't define this: shows not required to define all!
  // REQUIRE(a1 / a2 == A{1, 1});

  C c1{6};
  C c2{7};
  C c3{12};
  REQUIRE(c1 + 6 == C{12});
  REQUIRE(c1 - 3 == C{3});
  REQUIRE(6 + c1 == C{12});

  REQUIRE(c2 * 2 == C{14});
  REQUIRE(2 * c2 + 6 == C{20});
  REQUIRE(c2 / 2 == C{3});

  // REQUIRE(3 - c1 == C{-3});
}