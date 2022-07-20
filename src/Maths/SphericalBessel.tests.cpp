#include "SphericalBessel.hpp"
#include "catch2/catch.hpp"
#include "qip/Vector.hpp"
#include <iostream>
#include <vector>

TEST_CASE("Maths::SphericalBessel", "[jL][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "SphericalBessel\n";

  const auto r_list = qip::loglinear_range(1.0e-6, 200.0, 1.0, 75);

  for (int L = 0; L <= 15; ++L) {
    for (auto r : r_list) {
      REQUIRE(SphericalBessel::JL(L, r) ==
              Approx(SphericalBessel::exactGSL_JL(L, r)).margin(1.0e-9));
    }

    if (L < 4) {
      const auto jl = SphericalBessel::fillBesselVec(L, r_list);
      const auto jlkr = SphericalBessel::fillBesselVec_kr(L, 0.15, r_list);
      for (auto i = 0ul; i < r_list.size(); ++i) {
        const auto r = r_list.at(i);
        REQUIRE(jl.at(i) == Approx(SphericalBessel::JL(L, r)));
        REQUIRE(jlkr.at(i) == Approx(SphericalBessel::JL(L, 0.15 * r)));
      }
    }
  }
}