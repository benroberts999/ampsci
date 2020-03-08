#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

namespace Angular {

//******************************************************************************
// "Helper" functions
constexpr int twoj(int jindex) { return 2 * jindex + 1; }
constexpr int jindex(int twoj) { return (twoj - 1) / 2; }
constexpr int jindex_kappa(int ka) { return (ka > 0) ? ka - 1 : -ka - 1; }

constexpr int indexFromKappa(int ka) {
  return (ka < 0) ? -2 * ka - 2 : 2 * ka - 1;
}
constexpr int kappaFromIndex(int i) {
  return (i % 2 == 0) ? -(i + 2) / 2 : (i + 1) / 2;
}

template <typename T> inline T max4(T a, T b, T c, T d) {
  return std::max(std::max(a, b), std::max(c, d));
}
template <typename T> inline T min4(T a, T b, T c, T d) {
  return std::min(std::min(a, b), std::min(c, d));
}
inline int min_lambda_tj(int tja, int tjb, int tjc, int tjd) {
  return std::max(std::abs(tja - tjd), std::abs(tjb - tjc)) / 2;
}
inline int max_lambda_tj(int tja, int tjb, int tjc, int tjd) {
  return std::min((tja + tjd), (tjb + tjc)) / 2;
}

//******************************************************************************
class Ck_ab {
  // Lookup table for Ckab, and 3j symbols.
  // Ckab      = (-1)^{ja+1/2} * tildeCkab
  // tildeCkab = Sqrt([ja][jb]) * 3j(ja,jb,k, -1/2,1/2,0)pi

public:
  Ck_ab(const int in_max_K = 0, const int in_max_twoj = 0) {
    fill_maxK_twojmax(in_max_K, in_max_twoj);
  }

public:
  void fill_maxK_twojmax(const int in_max_K, const int in_max_twoj);

  // These will calculate values if they don't exist yet
  double get_tildeCkab_mutable(int k, int ka, int kb);
  double get_Ckab_mutable(int k, int ka, int kb);
  double get_3jkab_mutable(int k, int ka, int kb);

  // Use const versions if sure value already exists (will segfault otherwise)
  double get_tildeCkab(int k, int ka, int kb) const;
  double get_Ckab(int k, int ka, int kb) const;
  double get_3jkab(int k, int ka, int kb) const;

private:
  std::vector<std::vector<std::vector<double>>> m_3j_k_a_b = {};
  std::vector<std::vector<double>> m_Rjab_a_b = {}; // Sqrt([ja][jb])
  int m_max_jindex_sofar = -1;
  int m_max_k_sofar = -1;
};

//******************************************************************************

//******************************************************************************
class SixJTable_constk
//  { ja, jb, k}
//  { jc, jd, l}
// k is fixed
{
public:
  SixJTable_constk(int in_k, int tj_max = -1) : m_k(in_k), max_ji_sofar(-1) {
    fill(tj_max);
  }

private:
  const int m_k;
  int max_ji_sofar;
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>
      m_k_a_bcdl = {};

public:
  double get_6j(int tja, int tjb, int tjc, int tjd, int l) const;
  double get_6j(int tja, int tjb, int tjc, int tjd, int k, int l) const;
  //"safe" version that will allocate new 6j if not yet (NOT thresh-safe!)
  double get_6j_mutable(int tja, int tjb, int tjc, int tjd, int l);

  int get_k() const { return m_k; }
  int get_maxtj() const { return twoj(max_ji_sofar); }

public:
  void fill(const int tj_max);
};

//******************************************************************************
class SixJ
//  { ja, jb, k}
//  { jc, jd, l}
// any k
{
public:
  SixJ(int in_max_k = 0, int in_max_twoj = 1) { fill(in_max_k, in_max_twoj); }
  void fill(int in_max_k, int in_max_twoj);

  // Thread-safe, but will seg-fault if 6j doesn't exist
  double get_6j(int tja, int tjb, int tjc, int tjd, int k, int l) const;

  // Will calculate + store 6j if it doesn't exist, but not thread-safe
  double get_6j_mutable(int tja, int tjb, int tjc, int tjd, int k, int l);

private:
  int max_k_sofar = -1;
  int max_tj_sofar = -1;
  std::vector<SixJTable_constk> m_sixj_k = {};
};

} // namespace Angular

//
//
//******************************************************************************

//******************************************************************************

//******************************************************************************
//******************************************************************************
/*This was used to test method:*/

// double out = 0;
// {
//   ChronoTimer sw("table");
//   // auto sj = Angular::SixJ(2, 5);
//   for (int ttt = 0; ttt < 100; ttt++) {
//     for (int k = 0; k < 3; k++) {
//       for (int tja = 1; tja <= 5; tja += 2) {
//         for (int tjb = 1; tjb <= 5; tjb += 2) {
//           for (int tjc = 1; tjc <= 5; tjc += 2) {
//             for (int tjd = 1; tjd <= 5; tjd += 2) {
//               auto amd = std::abs(tja - tjd);
//               auto apd = tja + tjd;
//               auto bmc = std::abs(tjb - tjc);
//               auto bpc = tjb + tjc;
//               auto l_min = std::max(amd, bmc) / 2;
//               auto l_max = std::min(apd, bpc) / 2;
//               for (int l = l_min; l <= l_max; ++l) {
//                 out += Wigner::sixj_2(tja, tjb, 2 * k, tjc, tjd, 2 * l);
//                 // out += sj.get_6j(tja, tjb, tjc, tjd, k, l);
//                 //
//                 // std::cout << tja << "," << tjb << "," << tjc << "," << tjd
//                 // <<
//                 // "|"
//                 //           << k << "," << l << " :: ";
//                 // // auto v1 = sj.get_6j_mutable(tja, tjb, tjc, tjd, k, l);
//                 // auto v1 = sj.get_6j(tja, tjb, tjc, tjd, k, l);
//                 // auto v2 = Wigner::sixj_2(tja, tjb, 2 * k, tjc, tjd, 2 *
//                 l);
//                 // std::cout << v1 << " " << v2 << " = " << v1 - v2 << "\n";
//                 // if (std::abs(v1 - v2) > 1.e-10)
//                 //   std::cout << " ******************** \n ****************
//                 //   \n";
//               }
//             }
//           }
//         }
//       }
//     }
//   }
// }
// std::cout << out << "\n";
