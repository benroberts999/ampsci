#include "Wigner_369j.hpp"
#include <vector>

namespace Angular {

//******************************************************************************
// "Helper" functions
constexpr int twoj(int ji) { return 2 * ji + 1; }
constexpr int jindex(int twoj) { return (twoj - 1) / 2; }
inline int max4(int a, int b, int c, int d) {
  return std::max(std::max(a, b), std::max(c, d));
}
inline int min_lambda_tj(int tja, int tjb, int tjc, int tjd) {
  return std::max(abs(tja - tjd), abs(tjb - tjc)) / 2;
}
inline int max_lambda_tj(int tja, int tjb, int tjc, int tjd) {
  return std::min((tja + tjd), (tjb + tjc)) / 2;
}

//******************************************************************************
template <typename Real>
class SixJTable_constk
//  { ja, jb, k}
//  { jc, jd, l}
// k is fixed
{
public:
  SixJTable_constk(int in_k, int tj_max = -1) : k(in_k), max_ji_sofar(-1) {
    static_assert(std::is_floating_point<Real>::value,
                  "\nFAIL 25 in Angular: type of SixJTable must be "
                  "floating-point type!\n");
    fill(tj_max);
  }

private:
  const int k;
  int max_ji_sofar;
  std::vector<std::vector<std::vector<std::vector<std::vector<Real>>>>>
      m_k_a_bcdl;

public:
  Real value(int tja, int tjb, int tjc, int tjd, int l) {
    const auto lmin = min_lambda_tj(tja, tjb, tjc, tjd);
    if (l < lmin)
      return 0;
    if (l > max_lambda_tj(tja, tjb, tjc, tjd))
      return 0;

    auto a = jindex(tja);
    auto b = jindex(tjb);
    auto c = jindex(tjc);
    auto d = jindex(tjd);
    auto max = max4(a, b, c, d); //
    if (max == a) {
      return m_k_a_bcdl[a][b][c][d][l - lmin];
    } else if (max == b) {
      return m_k_a_bcdl[b][a][d][c][l - lmin];
    } else if (max == c) {
      return m_k_a_bcdl[c][d][a][b][l - lmin];
    } else {
      return m_k_a_bcdl[d][c][b][a][l - lmin];
    }
  }

  Real value(int tja, int tjb, int in_k, int tjc, int tjd, int l) {
    if (in_k == k)
      return value(tja, tjb, tjc, tjd, l);
    std::cerr << "\nFail 65 in Angular: Asked from from k; " << in_k << " " << k
              << "\n";
    std::abort();
    // return 0;
  }

public:
  void fill(const int tj_max) {

    const int ji_max = jindex(tj_max);
    if (ji_max <= max_ji_sofar)
      return;

    m_k_a_bcdl.reserve(ji_max + 1);
    for (int a = max_ji_sofar + 1; a <= ji_max; a++) {
      auto tja = twoj(a);
      std::vector<std::vector<std::vector<std::vector<Real>>>> ka_bcdl;
      ka_bcdl.reserve(a + 1);
      for (int b = 0; b <= a; b++) {
        auto tjb = twoj(b);

        std::vector<std::vector<std::vector<Real>>> ka_b_cdl;
        ka_b_cdl.reserve(a + 1);
        for (int c = 0; c <= a; c++) {
          auto tjc = twoj(c);
          const auto bmc = abs(tjb - tjc);
          const auto bpc = tjb + tjc;
          std::vector<std::vector<Real>> ka_bc_dl;
          ka_bc_dl.reserve(a + 1);
          for (int d = 0; d <= a; d++) {
            auto tjd = twoj(d);
            const auto amd = tja - tjd; // no abs, a>d
            const auto apd = tja + tjd;
            auto lambda_min = std::max(amd, bmc) / 2;
            auto lambda_max = std::min(apd, bpc) / 2;
            std::vector<Real> ka_bcd_l;
            ka_bcd_l.reserve(abs(lambda_max - lambda_min + 1));
            for (auto l = lambda_min; l <= lambda_max; l++) {
              ka_bcd_l.push_back(
                  Wigner::sixj_2(tja, tjb, 2 * k, tjc, tjd, 2 * l));
            }
            ka_bc_dl.push_back(ka_bcd_l);
          }
          ka_b_cdl.push_back(ka_bc_dl);
        }
        ka_bcdl.push_back(ka_b_cdl);
      }
      m_k_a_bcdl.push_back(ka_bcdl);
    }

    if (ji_max > max_ji_sofar)
      max_ji_sofar = ji_max;
  }
};

// //******************************************************************************
// template <typename Real> class SixJTable : public SixJTable_constk<Real> {
// public:
//   SixJTable(int k_max, int twoj_max = -1)
//       : // m_sixj.clear(); // shitty hack...
//         SixJTable_constk<Real>(0, -1) {
//     max_k = k_max; // xxx
//     for (int k = 0; k <= k_max; k++) {
//       m_sixj.push_back(SixJTable_constk<Real>(k, twoj_max));
//     }
//   }
//
// public: /*over-writes*/
//   Real value(int tja, int tjb, int tjc, int tjd, int l) = delete;
//   Real value(int tja, int tjb, int k, int tjc, int tjd, int l) {
//     if (k <= max_k)
//       return m_sixj[k].value(tja, tjb, tjc, tjd, l);
//     return 0;
//   }
//   void fill(const int tj_max) {
//     for (const auto &sixj : m_sixj) {
//       sixj.fill(tj_max);
//     }
//   }
//
// private:
//   int max_k;
//   std::vector<SixJTable_constk<Real>> m_sixj;
// }; // namespace Angular

} // namespace Angular
