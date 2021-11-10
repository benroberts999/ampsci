#pragma once
#include "BSpline.hpp"
#include "qip/Check.hpp"
#include <iostream>

bool BSplineClass(std::ostream &obuff) {
  bool pass = true;

  BSpline b(30, 4, 1.0e-4, 50.0);
  // b.print_knots();

  // Tests against mathematica, using the following MMA code
  // Knots are taken as output of BSpline class: b.print_knots()
  /*
  k = 4;
  d = k - 1;
  knots = {0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000100000000000, 0.000165649890835,
     0.000274398863336, 0.000454541417568, 0.000752947362000,
     0.001247256483196, 0.002066079002843, 0.003422457612770,
     0.005669297299420, 0.009391184787589, 0.015556487348734,
     0.025769304310905, 0.042686824459897, 0.070710678118655,
     0.117132161112086, 0.194029297014591, 0.321409318692165,
     0.532414185546389, 0.881943517146488, 1.460938473377593,
     2.420042986313319, 4.008798564982817, 6.640570446680343,
     11.000097695732780, 18.221649824700751, 30.184143042913092,
     50.000000000000000, 50.000000000000000, 50.000000000000000,
     50.000000000000000};
     BSplineBasis[{d, knots}, i, x]
     Evaluate[D[BSplineBasis[{d, knots}, i, x], {x, J}]]
     for J=(0,1,2)
     x = 0.01, 0.1, 1.0, 10.0
     i = [0, ..., N-1]
  */

  struct mma_data {
    double x;
    std::size_t i0;         // first non-zero spine index
    std::vector<double> b;  // non-zero splines, b(x)
    std::vector<double> b1; // non-zero spline derivatives, b'(x)
    std::vector<double> b2; // non-zero spline derivatives, b''(x)
  };
  // Data, from Mathematica
  std::vector<mma_data> test_data_list{
      {0.01,
       10,
       {0.231936, 0.658854, 0.109143, 0.0000671198},
       {-125.225, 52.1147, 72.7791, 0.33074},
       {45073.3, -70483.5, 24323.7, 1086.5}},
      {0.1,
       14,
       {0.0159261, 0.544434, 0.422133, 0.0175077},
       {-2.78881, -12.3797, 13.3753, 1.79325},
       {325.565, -319.263, -128.752, 122.451}},
      {1.0,
       19,
       {0.159858, 0.678553, 0.160998, 0.000590884},
       {-1.04043, 0.1041, 0.921316, 0.0150153},
       {4.51441, -6.81831, 2.04952, 0.254375}},
      {10.0,
       23,
       {0.00382507, 0.459167, 0.505112, 0.031896},
       {-0.0114741, -0.145493, 0.128484, 0.0284834},
       {0.022946, -0.00870835, -0.0311949, 0.0169573}}};

  // Compare splines to MMA, for 0th, 1st, and 2nd order derivatives:
  double max_eps = 0.0;
  std::size_t max_deli = 0;
  for (const auto &[x, ti0, tb, tb1, tb2] : test_data_list) {
    // evaulate the splines:
    const auto [i0, bij] = b.get_nonzero(x, 2);
    for (std::size_t i = 0; i < b.K(); ++i) {
      const auto deli0 = i0 > ti0 ? i0 - ti0 : ti0 - i0;
      const auto del0 = std::abs((bij[i][0] - tb[i]) / tb[i]);
      const auto del1 = std::abs((bij[i][1] - tb1[i]) / tb1[i]);
      const auto del2 = std::abs((bij[i][2] - tb2[i]) / tb2[i]);
      const auto eps = std::max({del0, del1, del2});
      printf("%2i/%2i %9.3e [%9.3e], %10.3e [%10.3e], %10.3e [%10.3e]\n",
             int(i0), int(ti0), bij[i][0], tb[i], bij[i][1], tb1[i], bij[i][2],
             tb2[i]);
      if (eps > max_eps) {
        max_eps = eps;
      }
      if (deli0 > max_deli) {
        max_deli = deli0;
      }
    }
  }

  pass &= qip::check_value(&obuff, "Spline i0", int(max_deli), 0, 0);
  pass &= qip::check_value(&obuff, "Spline values", max_eps, 0.0, 1.0e-5);

  //****************************************************************************
  { // Check that the splines sum to 1
    // This time for N=90, K=9
    BSpline b2(90, 9, 1.0e-5, 70.0);
    double delta = 0.0;
    for (auto &x : {1.0e-5, 0.01, 0.1, 10.0, 25.0, 70.0}) {
      const auto [i0, bij] = b2.get_nonzero(x, 0);
      double sum = 0.0;
      for (auto i = 0ul; i < b2.K(); ++i) {
        sum += bij[i][0];
      }
      delta += std::abs(1.0 - sum);
    }
    pass &= qip::check_value(&obuff, "Spline Sum", delta, 0.0, 1.0e-14);
  }
  { // Check that the splines sum to 1
    // This time for N=90, K=9
    BSpline b2(90, 9, 1.0e-5, 70.0, BSpline::KnotDistro::loglinear);
    double delta = 0.0;
    for (auto &x : {1.0e-5, 0.01, 0.1, 10.0, 25.0, 70.0}) {
      const auto [i0, bij] = b2.get_nonzero(x, 0);
      double sum = 0.0;
      for (auto i = 0ul; i < b2.K(); ++i) {
        sum += bij[i][0];
      }
      delta += std::abs(1.0 - sum);
    }
    pass &= qip::check_value(&obuff, "Spline Sum LL", delta, 0.0, 1.0e-14);
  }

  return pass;
}
