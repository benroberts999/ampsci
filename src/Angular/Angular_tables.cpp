#include "Angular/Angular_tables.hpp"
#include "Angular/Angular_369j.hpp"
#include <cassert>
// NB: These headers MUST be in this order; otherwise fails to compile
// on macOS... what :\ At least it works now....
// Perhaps something to do with cmath vs math.h conflict???
// math.h might be included by gsl, which is in Angular_369j...
// Can only reproduce issue on mac..yay

#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace Angular {
//******************************************************************************
void Ck_ab::fill_maxK_twojmax(const int in_max_K, const int in_max_twoj) {

  auto max_jindex = (in_max_twoj - 1) / 2;
  if (max_jindex < m_max_jindex_sofar)
    max_jindex = m_max_jindex_sofar;
  auto max_k = in_max_K > m_max_k_sofar ? in_max_K : m_max_k_sofar;

  m_Rjab_a_b.resize(max_jindex + 1);
  for (int jia = m_max_jindex_sofar + 1; jia <= max_jindex; jia++) {
    m_Rjab_a_b[jia].reserve(jia + 1);
    for (int jib = 0; jib <= jia; jib++) {
      auto tjpqtjpq = static_cast<double>((2 * jia + 2) * (2 * jib + 2));
      m_Rjab_a_b[jia].push_back(std::sqrt(tjpqtjpq));
    }
  }

  m_3j_k_a_b.resize(max_k + 1);
  for (int k = 0; k <= max_k; k++) {
    m_3j_k_a_b[k].resize(max_jindex + 1);
    auto beg_jia = (k <= m_max_k_sofar) ? m_max_jindex_sofar + 1 : 0;
    for (int jia = beg_jia; jia <= max_jindex; jia++) {
      m_3j_k_a_b[k][jia].reserve(jia + 1);
      for (int jib = 0; jib <= jia; jib++) {
        auto tja = 2 * jia + 1;
        auto tjb = 2 * jib + 1;
        m_3j_k_a_b[k][jia].push_back(
            Angular::special_threej_2(tja, tjb, 2 * k));
      }
    }
  }

  m_max_jindex_sofar = max_jindex;
  m_max_k_sofar = in_max_K;
}

//******************************************************************************
double Ck_ab::get_tildeCkab_mutable(int k, int ka, int kb) {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  // parity:
  auto pi_ok = Angular::evenQ(Angular::l_k(ka) + Angular::l_k(kb) + k);
  if (!pi_ok)
    return 0;

  auto maxji = std::max(jia, jib);
  if (maxji > m_max_jindex_sofar || k > m_max_k_sofar)
    fill_maxK_twojmax(k, 2 * maxji + 1);
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] * m_Rjab_a_b[jia][jib]
                     : m_3j_k_a_b[k][jib][jia] * m_Rjab_a_b[jib][jia];
}

double Ck_ab::get_tildeCkab(int k, int ka, int kb) const {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  // parity:
  auto pi_ok = Angular::evenQ(Angular::l_k(ka) + Angular::l_k(kb) + k);
  if (!pi_ok)
    return 0;
  assert(std::max(jia, jib) <= m_max_jindex_sofar);
  assert(k <= m_max_k_sofar);
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] * m_Rjab_a_b[jia][jib]
                     : m_3j_k_a_b[k][jib][jia] * m_Rjab_a_b[jib][jia];
}

//******************************************************************************
double Ck_ab::get_Ckab_mutable(int k, int ka, int kb) {
  auto s = Angular::evenQ_2(Angular::twoj_k(ka) + 1) ? 1.0 : -1.0;
  return s * get_tildeCkab_mutable(k, ka, kb);
}

double Ck_ab::get_Ckab(int k, int ka, int kb) const {
  auto s = Angular::evenQ_2(Angular::twoj_k(ka) + 1) ? 1.0 : -1.0;
  return s * get_tildeCkab(k, ka, kb);
}

//******************************************************************************
double Ck_ab::get_3jkab_mutable(int k, int ka, int kb) {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  auto maxji = std::max(jia, jib);
  if (maxji > m_max_jindex_sofar || k > m_max_k_sofar)
    fill_maxK_twojmax(k, 2 * maxji + 1);
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] : m_3j_k_a_b[k][jib][jia];
}

double Ck_ab::get_3jkab(int k, int ka, int kb) const {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  assert(std::max(jia, jib) <= m_max_jindex_sofar);
  assert(k <= m_max_k_sofar);
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] : m_3j_k_a_b[k][jib][jia];
}

double Ck_ab::get_Lambdakab(int k, int ka, int kb) const {
  const auto pi_ok = Angular::evenQ(Angular::l_k(ka) + Angular::l_k(kb) + k);
  if (!pi_ok)
    return 0.0;
  const auto tjs = get_3jkab(k, ka, kb);
  return tjs * tjs;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
double SixJTable_constk::get_6j(int tja, int tjb, int tjc, int tjd,
                                int l) const {
  const auto lmin = min_lambda_tj(tja, tjb, tjc, tjd);
  if (l < lmin)
    return 0.0;
  if (l > max_lambda_tj(tja, tjb, tjc, tjd))
    return 0.0;

  const auto a = jindex(tja);
  const auto b = jindex(tjb);
  const auto c = jindex(tjc);
  const auto d = jindex(tjd);
  const auto max = max4(a, b, c, d); //
  assert(max <= max_ji_sofar);

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

//******************************************************************************
double SixJTable_constk::get_6j(int tja, int tjb, int tjc, int tjd, int in_k,
                                int l) const {
  if (in_k == m_k)
    return get_6j(tja, tjb, tjc, tjd, l);
  std::cerr << "\nFail 65 in Angular: Asked from from k; " << in_k << " " << m_k
            << "\n";
  std::abort();
}

//******************************************************************************
double SixJTable_constk::get_6j_mutable(int tja, int tjb, int tjc, int tjd,
                                        int l) {
  auto maxtj = max4(tja, tjb, tjc, tjd);
  if (jindex(maxtj) > max_ji_sofar)
    fill(maxtj);
  return get_6j(tja, tjb, tjc, tjd, l);
}

//******************************************************************************
void SixJTable_constk::fill(const int tj_max) {

  const int ji_max = jindex(tj_max);
  if (ji_max <= max_ji_sofar)
    return;

  m_k_a_bcdl.reserve(ji_max + 1);
  for (int a = max_ji_sofar + 1; a <= ji_max; a++) {
    auto tja = twoj(a);
    std::vector<std::vector<std::vector<std::vector<double>>>> ka_bcdl;
    ka_bcdl.reserve(a + 1);
    for (int b = 0; b <= a; b++) {
      auto tjb = twoj(b);

      std::vector<std::vector<std::vector<double>>> ka_b_cdl;
      ka_b_cdl.reserve(a + 1);
      for (int c = 0; c <= a; c++) {
        auto tjc = twoj(c);
        const auto bmc = std::abs(tjb - tjc);
        const auto bpc = tjb + tjc;
        std::vector<std::vector<double>> ka_bc_dl;
        ka_bc_dl.reserve(a + 1);
        for (int d = 0; d <= a; d++) {
          auto tjd = twoj(d);
          const auto amd = tja - tjd; // no abs, a>d
          const auto apd = tja + tjd;
          auto lambda_min = std::max(amd, bmc) / 2;
          auto lambda_max = std::min(apd, bpc) / 2;
          std::vector<double> ka_bcd_l;
          ka_bcd_l.reserve(std::abs(lambda_max - lambda_min + 1));
          for (auto l = lambda_min; l <= lambda_max; l++) {
            ka_bcd_l.push_back(
                Angular::sixj_2(tja, tjb, 2 * m_k, tjc, tjd, 2 * l));
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

//******************************************************************************
//******************************************************************************
void SixJ::fill(int in_max_k, int in_max_twoj) {
  // new max_tj for each existing 6js:
  //(each sixj_constk keeps track of max_2j)
  if (in_max_twoj < max_tj_sofar)
    in_max_twoj = max_tj_sofar;
  for (auto &sixjk : m_sixj_k) {
    sixjk.fill(in_max_twoj);
  }
  // Now, fill new k's
  m_sixj_k.reserve(in_max_k + 1);
  for (int k = max_k_sofar + 1; k <= in_max_k; ++k) {
    m_sixj_k.emplace_back(k, in_max_twoj);
  }
  max_tj_sofar = in_max_twoj;
  max_k_sofar = in_max_k;
}

//******************************************************************************
double SixJ::get_6j(int tja, int tjb, int tjc, int tjd, int k, int l) const {
  assert(k <= max_k_sofar);
  return m_sixj_k[k].get_6j(tja, tjb, tjc, tjd, l);
}

//******************************************************************************
double SixJ::get_6j_mutable(int tja, int tjb, int tjc, int tjd, int k, int l) {
  if (k > max_k_sofar) {
    fill(k, max4(tja, tjb, tjc, tjd));
  }
  return m_sixj_k[k].get_6j_mutable(tja, tjb, tjc, tjd, l);
}

} // namespace Angular
