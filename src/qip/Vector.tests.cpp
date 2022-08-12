#include "Vector.hpp"
#include "catch2/catch.hpp"
#include <cassert>
#include <iostream>

TEST_CASE("qip::Vector", "[qip][Vector][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "qip::Vector\n";

  const std::vector a{1, 2, 3, 4, 5, 6};
  const std::vector b{1, 2, 4, 4, 5, 6};

  const auto [del, it] = qip::compare(a, b);
  REQUIRE(*it == 3);
  REQUIRE(del == -1);

  const auto chi2 = [](auto x, auto y) { return (x - y) * (x - y); };
  const auto [del2, it2] = qip::compare(a, b, chi2);
  REQUIRE(*it2 == 3);
  REQUIRE(del2 == 1);

  const std::vector x{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  const std::vector y{1.0, 2.0, 4.0, 4.0, 5.0, 6.0};

  const auto [eps3, it3] = qip::compare_eps(x, y);
  REQUIRE(*it3 == 3.0);
  REQUIRE(eps3 == Approx(-1.0 / 4.0));

  const std::vector c{2, 4, 6, 8};
  const auto d = qip::add(a, b, c);
  REQUIRE(d == std::vector{4, 8, 13, 16, 10, 12});

  std::vector e{1, 1};
  qip::add(&e, d);
  REQUIRE(e == std::vector{5, 9, 13, 16, 10, 12});

  const auto d2 = qip::multiply(a, b, c);
  REQUIRE(d2 == std::vector{2, 16, 72, 128, 0, 0});

  std::vector e2{2, 2};
  qip::multiply(&e2, d2);
  REQUIRE(e2 == std::vector{4, 32, 0, 0, 0, 0});

  const auto d3 = qip::compose(chi2, a, b, c);
  REQUIRE(d3 == std::vector{4, 16, 25, 64, 0, 0});

  std::vector e3{2, 2};
  qip::compose(chi2, &e3, d3);
  REQUIRE(e3 == std::vector{4, 196, 625, 4096, 0, 0});

  const auto d4 = qip::scale(d3, 2);
  REQUIRE(d4 == std::vector{8, 32, 50, 128, 0, 0});

  std::vector e4{2, 2};
  qip::scale(&e4, 2);
  REQUIRE(e4 == std::vector{4, 4});

  const auto f = qip::inner_product(a, b, d);
  REQUIRE(f == (4 + 32 + 156 + 256 + 250 + 432));

  const auto f2 = qip::inner_product_sub(1, 5, a, b, d);
  REQUIRE(f2 == (32 + 156 + 256 + 250));

  const auto f3 = qip::inner_product(x, y);
  REQUIRE(f3 == Approx(1.0 + 4.0 + 12.0 + 16.0 + 25.0 + 36.0));

  auto square = [](auto t) { return t * t; };
  auto g = qip::apply_to(square, d4);
  REQUIRE(g == std::vector{64, 1024, 2500, 16384, 0, 0});

  //============================================================================
  {
    using namespace qip::overloads;
    //
    const auto d0 = a + b + c;
    REQUIRE(d0 == d);
    auto a0 = a;
    a0 += (b + c);
    REQUIRE(a0 == d);
    a0 -= (b + c);
    REQUIRE(a0 == a);
    const auto a00 = d0 - b - c;
    REQUIRE(a00 == a);

    const auto dx0 = a * 2;
    REQUIRE(dx0 == std::vector{2, 4, 6, 8, 10, 12});
    auto ax0 = a;
    ax0 *= 2;
    REQUIRE(ax0 == dx0);
    ax0 /= 2;
    REQUIRE(ax0 == a);
    const auto ax00 = dx0 / 2;
    REQUIRE(ax00 == a);

    const std::vector w{1, 2};
    const std::vector z{1, 2, 3, 4, 5, 6};

    REQUIRE(w + z == std::vector{2, 4, 3, 4, 5, 6});
    REQUIRE(z + w == std::vector{2, 4, 3, 4, 5, 6});
    auto wz = w;
    wz += z;
    REQUIRE(wz == std::vector{2, 4, 3, 4, 5, 6});
  }

  //============================================================================

  const auto ur = qip::uniform_range(1, 10, 10);
  REQUIRE(ur == std::vector{1, 2, 3, 4, 5, 6, 7, 8, 9, 10});

  const auto ur2 = qip::uniform_range(2, 10, 5);
  REQUIRE(ur2 == std::vector{2, 4, 6, 8, 10});

  const auto ur3 = qip::uniform_range(2.0, 10.0, 5);
  REQUIRE(ur3.front() == 2.0);
  REQUIRE(ur3.back() == 10.0);

  const auto lr = qip::logarithmic_range(1, 10000, 5);
  REQUIRE(lr == std::vector{1, 10, 100, 1000, 10000});

  const auto lr2 = qip::logarithmic_range(1.0, 10000.0, 5);
  REQUIRE(lr2.front() == 1.0);
  REQUIRE(lr2.back() == 10000.0);

  const auto llr = qip::loglinear_range(1.0, 100.0, 10.0, 50);
  REQUIRE(llr.front() == 1.0);
  REQUIRE(llr.back() == 100.0);
}