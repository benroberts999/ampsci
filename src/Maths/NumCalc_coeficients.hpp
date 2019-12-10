#pragma once
#include <array>

namespace NumCalc {
//******************************************************************************
// Define the coeficients for quadrature integration:
template <std::size_t N> struct QintCoefs {};

template <> struct QintCoefs<13> {
  static constexpr std::size_t N = 13;
  static constexpr std::array<double, N> cq{
      {1382741929621, 9535909891802, -5605325192308, 28323664941310,
       -32865015189975, 53315213499588, -41078125154304, 39022895874876,
       -13155015007785, 12465244770050, 3283609164916, 5551687979302,
       5206230892907}};
  static constexpr double dq_inv = 1.0 / 5230697472000;
};
template <> struct QintCoefs<11> {
  static constexpr std::size_t N = 11;
  static constexpr std::array<double, N> cq{
      {262747265, 1637546484, -454944189, 3373884696, -2145575886, 3897945600,
       -1065220914, 1942518504, 636547389, 1021256716, 952327935}};
  static constexpr double dq_inv = 1.0 / 958003200;
};
template <> struct QintCoefs<9> {
  static constexpr std::size_t N = 9;
  static constexpr std::array<double, N> cq{{2082753, 11532470, 261166,
                                             16263486, -1020160, 12489922,
                                             5095890, 7783754, 7200319}};
  static constexpr double dq_inv = 1.0 / 7257600;
};
template <> struct QintCoefs<7> {
  static constexpr std::size_t N = 7;
  static constexpr std::array<double, N> cq{
      {36799, 176648, 54851, 177984, 89437, 130936, 119585}};
  static constexpr double dq_inv = 1.0 / 120960;
};
template <> struct QintCoefs<5> {
  static constexpr std::size_t N = 5;
  static constexpr std::array<double, N> cq{{475, 1902, 1104, 1586, 1413}};
  static constexpr double dq_inv = 1.0 / 1440;
};
template <> struct QintCoefs<3> {
  static constexpr std::size_t N = 3;
  static constexpr std::array<double, N> cq{{9, 28, 23}};
  static constexpr double dq_inv = 1.0 / 24;
};
template <> struct QintCoefs<1> {
  static constexpr std::size_t N = 1;
  static constexpr std::array<double, N> cq{{1}};
  static constexpr double dq_inv = 1.0 / 2;
};

} // namespace NumCalc
