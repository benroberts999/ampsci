#include "Methods.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <vector>

TEST_CASE("qip::Methods", "[qip][Methods][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "qip::Methods\n";

  {
    const auto y_x2 = [](auto x) { return x * x; };
    const auto y_x2_root = 0.0;

    const auto [x0, dx0] = qip::Newtons(y_x2, 0.01);
    const auto error0 = std::abs(x0 - y_x2_root);
    REQUIRE(x0 == Approx(y_x2_root).margin(2.0 * 1.0e-6));
    REQUIRE(error0 <= 2.0 * dx0);

    const double delta_target = 1.0e-9;
    const auto [x1, dx1] = qip::Newtons(y_x2, 0.01, delta_target);
    const auto error1 = std::abs(x1 - y_x2_root);

    REQUIRE(x1 == Approx(0.0).margin(2.0 * delta_target));
    REQUIRE(error1 <= 2.0 * dx1);

    {
      const auto [dydx, err] = qip::derivative(y_x2, 0.0, delta_target);
      const auto err_actual = std::abs(dydx - 0.0);
      REQUIRE(dydx == Approx(0.0).margin(2.0 * delta_target));
      REQUIRE(err_actual <= 2.0 * err);
    }
    {
      const auto [dydx, err] = qip::derivative(y_x2, 12.0, delta_target);
      const auto err_actual = std::abs(dydx - 24.0);
      REQUIRE(dydx == Approx(24.0).margin(2.0 * delta_target));
      REQUIRE(err_actual <= 2.0 * err);
    }
  }

  {
    const auto y_x3 = [](auto x) { return (x * x * (x - 10.5)); };
    const auto y_x3_root = 10.5;

    const double delta_target = 1.0e-9;
    const double initial_guess = 8.0;
    const auto [x0, dx0] = qip::Newtons(y_x3, initial_guess, delta_target);
    const auto error = std::abs(x0 - y_x3_root);
    REQUIRE(x0 == Approx(y_x3_root).margin(2.0 * delta_target));
    REQUIRE(error < 2.0 * dx0);

    {
      const auto [dydx, err] = qip::derivative(y_x3, 0.0, delta_target);
      const auto err_actual = std::abs(dydx - 0.0);
      REQUIRE(dydx == Approx(0.0).margin(2.0 * delta_target));
      REQUIRE(err_actual <= 2.0 * err);
    }
    {
      const auto [dydx, err] = qip::derivative(y_x3, 3.5, delta_target);
      const auto err_actual = std::abs(dydx - -36.75);
      REQUIRE(dydx == Approx(-36.75).margin(2.0 * delta_target));
      REQUIRE(err_actual <= 2.0 * err);
    }
  }
}