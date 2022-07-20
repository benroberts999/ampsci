#include "Grid.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <vector>

TEST_CASE("Maths::Grid", "[Grid][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Grid\n";

  for (const auto &gr_type :
       {GridType::loglinear, GridType::logarithmic, GridType::linear}) {

    auto gt_str = GridParameters::parseType(gr_type);
    REQUIRE(GridParameters::parseType(gt_str) == gr_type);

    const auto r0 = 1.0e-6;
    const auto rmax = 100.0;
    const auto num_pts = 1000ul;
    const auto b = 4.0;

    // test the various constructors
    const Grid g1(r0, rmax, num_pts, gr_type, b);
    Grid g2(GridParameters{num_pts, r0, rmax, b, gr_type});
    const Grid g3(GridParameters{num_pts, r0, rmax, b, gt_str});
    Grid g4(GridParameters{0, r0, rmax, b, gr_type, g1.du()});
    const auto g5 = g1;
    Grid g6(1.0e-2, 66.0, 123, GridType::linear, 0.0);
    g6 = g1;
    const Grid g7(g1.params());

    REQUIRE(g1.du() == Approx(g4.du()));
    REQUIRE(g1.du() == Approx(g7.du()));

    REQUIRE(g1.r() == g2.r());
    REQUIRE(g1.r() == g3.r());
    REQUIRE(g1.r() == g4.r());
    REQUIRE(g1.r() == g5.r());
    REQUIRE(g1.r() == g6.r());

    REQUIRE(g1.drdu() == g2.drdu());
    REQUIRE(g1.drdu() == g3.drdu());
    REQUIRE(g1.drduor() == g2.drduor());
    REQUIRE(g1.drduor() == g3.drduor());
    REQUIRE(g1.drduor() == g5.drduor());
    REQUIRE(g1.drduor() == g6.drduor());

    REQUIRE(g1.r(3) == g2.r(3));
    REQUIRE(g1.r(3) == g3.r(3));
    REQUIRE(g1.r(3) == g4.r(3));
    REQUIRE(g1.drdu(3) == g2.drdu(3));
    REQUIRE(g1.drdu(3) == g3.drdu(3));
    REQUIRE(g1.drduor(3) == g2.drduor(3));
    REQUIRE(g1.drduor(3) == g3.drduor(3));

    REQUIRE(g1.r0() == Approx(r0));
    REQUIRE(g1.r0() == g2.r0());
    REQUIRE(g1.r0() == g3.r0());
    REQUIRE(g1.r0() == g4.r0());

    REQUIRE(g1.rmax() == Approx(rmax));
    REQUIRE(g1.rmax() == g2.rmax());
    REQUIRE(g1.rmax() == g3.rmax());
    REQUIRE(g1.rmax() == g4.rmax());

    REQUIRE(g1.num_points() == num_pts);
    REQUIRE(g1.num_points() == g2.num_points());
    REQUIRE(g1.num_points() == g3.num_points());
    REQUIRE(g1.num_points() == g4.num_points());
    REQUIRE(g1.num_points() == g7.num_points());

    REQUIRE(g1.type() == gr_type);
    REQUIRE(g1.type() == g2.type());
    REQUIRE(g1.type() == g3.type());
    REQUIRE(g1.type() == g4.type());
    REQUIRE(g1.type() == g7.type());

    const auto t_b = gr_type == GridType::loglinear ? b : 0.0;
    REQUIRE(g1.loglin_b() == Approx(t_b));
    REQUIRE(g1.loglin_b() == g2.loglin_b());
    REQUIRE(g1.loglin_b() == g3.loglin_b());
    REQUIRE(g1.loglin_b() == g4.loglin_b());
    REQUIRE(g1.loglin_b() == g7.loglin_b());

    REQUIRE(g1.gridParameters() == g2.gridParameters());

    auto i01 = g1.getIndex(r0, false);
    auto i02 = g1.getIndex(r0, true);
    REQUIRE(i01 == 0);
    REQUIRE(i02 == 0);
    auto if1 = g1.getIndex(rmax, false);
    auto if2 = g1.getIndex(rmax, true);
    REQUIRE(if1 == 999);
    REQUIRE(if2 == 999);

    auto ix0 = g1.getIndex(g1.r(60), false);
    auto ix1 = g1.getIndex(g1.r(60), true);
    REQUIRE(ix0 == Approx(60).margin(1));
    REQUIRE(ix1 == 60);

    const auto delta = 0.2 * (g1.r(61) - g1.r(60));

    auto ix3 = g1.getIndex(g1.r(60) + delta, false);
    auto ix4 = g1.getIndex(g1.r(60) + delta, true);
    REQUIRE(ix3 == Approx(60).margin(1));
    REQUIRE(ix4 == 60);

    auto ix5 = g1.getIndex(g1.r(60) - delta, false);
    auto ix6 = g1.getIndex(g1.r(60) - delta, true);
    REQUIRE(ix5 == Approx(60).margin(1));
    REQUIRE(ix6 == 60);
  }
}