#pragma once
#include <iostream>
#include <vector>

namespace NumCalc {

//******************************************************************************
inline std::vector<double> derivative(const std::vector<double> &f,
                                      const std::vector<double> &drdt,
                                      const double dt, const int order = 1)
/*
df/dr = df/dt * dt/dr = (df/dt) / (dr/dt) = (df/di) * (di/dt) / (dr/dt) =
= (df/di)  / (dt * dr/dt)
coeficients from: http://en.wikipedia.org/wiki/Finite_difference_coefficient
*/
{

  size_t ngp = f.size();
  std::vector<double> df(ngp);

  df[0] = (f[1] - f[0]) / (dt * drdt[0]);
  df[ngp - 1] = (f[ngp - 1] - f[ngp - 2]) / (dt * drdt[ngp - 1]);

  df[1] = (f[2] - f[0]) / (2 * dt * drdt[1]);
  df[ngp - 2] = (f[ngp - 1] - f[ngp - 3]) / (2 * dt * drdt[ngp - 2]);

  df[2] = (f[0] - 8 * f[1] + 8 * f[3] - f[4]) / (12 * dt * drdt[2]);
  df[ngp - 3] = (f[ngp - 5] - 8 * f[ngp - 4] + 8 * f[ngp - 2] - f[ngp - 1]) /
                (12 * dt * drdt[ngp - 3]);

  df[3] = (-1 * f[0] + 9 * f[1] - 45 * f[2] + 45 * f[4] - 9 * f[5] + 1 * f[6]) /
          (60 * dt * drdt[3]);
  df[ngp - 4] = (-1 * f[ngp - 7] + 9 * f[ngp - 6] - 45 * f[ngp - 5] +
                 45 * f[ngp - 3] - 9 * f[ngp - 2] + 1 * f[ngp - 1]) /
                (60 * dt * drdt[ngp - 4]);

  for (size_t i = 4; i < (ngp - 4); i++) {
    df[i] = ((1. / 8) * f[i - 4] - (4. / 3) * f[i - 3] + 7 * f[i - 2] -
             28 * f[i - 1] - (1. / 8) * f[i + 4] + (4. / 3) * f[i + 3] -
             7 * f[i + 2] + 28 * f[i + 1]) /
            (35 * dt * drdt[i]);
  }

  if (order > 1)
    df = derivative(df, drdt, dt, order - 1);

  return df;
}

//******************************************************************************
// static size_t Nquad = 9;
// static double cq[9] = {2082753.,  11532470., 261166.,  16263486., -1020160.,
//                        12489922., 5095890.,  7783754., 7200319.};
// static double dq_inv = 1. / 7257600.;
static size_t Nquad = 3;
static double cq[3] = {9., 28., 23.};
static double dq_inv = 1. / 24.;
//******************************************************************************
inline double integrate(const std::vector<double> &f1, const double dt = 1.,
                        size_t beg = 0, size_t end = 0)
/*
Note: includes no safety checks!
Integrates from (point) beg to end-1 (i.e., not including end)
Require:
  * (beg-end) > 2*Nquad
  * end - Nquad > Nquad
  * beg + 2*Nquad
  Really, we probably need : (beg-end) >> 2*Nquad
*/
{

  if (end == 0)
    end = f1.size();

  // if (end - beg < 2 * Nquad)
  //   std::cerr << "\nFAIL 71 in INT: interval too small\n";

  double Rint_s = 0;
  for (size_t i = 0; i < Nquad; i++)
    Rint_s += cq[i] * f1[beg + i];

  double Rint_m = 0;
  for (auto i = beg + Nquad; i < end - Nquad; i++)
    Rint_m += f1[i];

  double Rint_e = 0;
  for (size_t i = 0; i < Nquad; i++)
    Rint_e += cq[i] * f1[end - i - 1];

  return (Rint_m + dq_inv * (Rint_s + Rint_e)) * dt;

} // END integrate1

//******************************************************************************
inline double integrate(const std::vector<double> &f1,
                        const std::vector<double> &f2, const double dt = 1.,
                        size_t beg = 0, size_t end = 0)
// Copy-paste from above
{

  if (end == 0)
    end = f1.size();

  // if (end - beg < 2 * Nquad)
  //   std::cerr << "\nFAIL 71 in INT: interval too small\n";

  double Rint_s = 0;
  for (size_t i = 0; i < Nquad; i++)
    Rint_s += cq[i] * f1[beg + i] * f2[beg + i];

  double Rint_m = 0;
  for (auto i = beg + Nquad; i < end - Nquad; i++)
    Rint_m += f1[i] * f2[i];

  double Rint_e = 0;
  for (size_t i = 0; i < Nquad; i++)
    Rint_e += cq[i] * f1[end - i - 1] * f2[end - i - 1];

  return (Rint_m + dq_inv * (Rint_s + Rint_e)) * dt;

} // END integrate 2

//******************************************************************************
inline double integrate(const std::vector<double> &f1,
                        const std::vector<double> &f2,
                        const std::vector<double> &f3, const double dt = 1.,
                        size_t beg = 0, size_t end = 0)
// Copy-paste from above
{

  if (end == 0)
    end = f1.size();

  // if (end - beg < 2 * Nquad)
  //   std::cerr << "\nFAIL 71 in INT: interval too small\n";

  double Rint_s = 0;
  for (size_t i = 0; i < Nquad; i++)
    Rint_s += cq[i] * f1[beg + i] * f2[beg + i] * f3[beg + i];

  double Rint_m = 0;
  for (auto i = beg + Nquad; i < end - Nquad; i++)
    Rint_m += f1[i] * f2[i] * f3[i];

  double Rint_e = 0;
  for (size_t i = 0; i < Nquad; i++)
    Rint_e += cq[i] * f1[end - i - 1] * f2[end - i - 1] * f3[end - i - 1];

  return (Rint_m + dq_inv * (Rint_s + Rint_e)) * dt;

} // END integrate 3

//******************************************************************************
inline double
integrate(const std::vector<double> &f1, const std::vector<double> &f2,
          const std::vector<double> &f3, const std::vector<double> &f4,
          const double dt = 1., size_t beg = 0, size_t end = 0)
// Copy-paste from above
{

  if (end == 0)
    end = f1.size();

  // if (end - beg < 2 * Nquad)
  //   std::cerr << "\nFAIL 71 in INT: interval too small\n";

  double Rint_s = 0;
  for (size_t i = 0; i < Nquad; i++)
    Rint_s += cq[i] * f1[beg + i] * f2[beg + i] * f3[beg + i] * f4[beg + i];

  double Rint_m = 0;
  for (auto i = beg + Nquad; i < end - Nquad; i++)
    Rint_m += f1[i] * f2[i] * f3[i] * f4[i];

  double Rint_e = 0;
  for (size_t i = 0; i < Nquad; i++)
    Rint_e += cq[i] * f1[end - i - 1] * f2[end - i - 1] * f3[end - i - 1] *
              f4[end - i - 1];

  return (Rint_m + dq_inv * (Rint_s + Rint_e)) * dt;

} // END integrate 4

} // namespace NumCalc

//******************************************************************************

//******************************************************************************

//******************************************************************************

/*
if (nquad == 1) {
  double c[1] = {1};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 2) {
  double c[2] = {1, 2};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 3) {
  double c[3] = {9, 28, 23};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 4) {
  double c[4] = {8, 31, 20, 25};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 5) {
  double c[5] = {475, 1902, 1104, 1586, 1413};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 6) {
  double c[6] = {459, 1982, 944, 1746, 1333, 1456};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 7) {
  double c[7] = {36799, 176648, 54851, 177984, 89437, 130936, 119585};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 8) {
  double c[8] = {35584, 185153, 29336, 220509, 46912, 156451, 111080, 122175};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 9) {
  double c[9] = {2082753,  11532470, 261166,  16263486, -1020160,
                 12489922, 5095890,  7783754, 7200319};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 10) {
  double c[10] = {2034625,  11965622, -1471442, 20306238, -7084288,
                  18554050, 1053138,  9516362,  6767167,  7305728};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 11) {
  double c[11] = {262747265,   1637546484, -454944189,  3373884696,
                  -2145575886, 3897945600, -1065220914, 1942518504,
                  636547389,   1021256716, 952327935};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 12) {
  double c[12] = {257696640,   1693103359, -732728564,  4207237821,
                  -3812282136, 6231334350, -3398609664, 3609224754,
                  -196805736,  1299041091, 896771060,   963053825};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 13) {
  double c[13] = {1382741929621,   9535909891802,   -5605325192308,
                  28323664941310,  -32865015189975, 53315213499588,
                  -41078125154304, 39022895874876,  -13155015007785,
                  12465244770050,  3283609164916,   5551687979302,
                  5206230892907};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else if (nquad == 14) {
  double c[14] = {1360737653653,   9821965479386,   -7321658717812,
                  34616887868158,  -48598072507095, 81634716670404,
                  -78837462715392, 76782233435964,  -41474518178601,
                  28198302087170,  -3009613761932,  7268021504806,
                  4920175305323,   5252701747968};
  for (int i = 0; i < nquad; i++) {
    ic[i] = c[i];
  }
} else {
  // printf("FAILURE: Wrong order for integration.. check nquad\n");
  return 1; // XXX
}

double dd[14] = {2,
                 2,
                 24,
                 24,
                 1440,
                 1440,
                 120960,
                 120960,
                 7257600,
                 7257600,
                 958003200,
                 958003200,
                 5230697472000,
                 5230697472000};

*/
