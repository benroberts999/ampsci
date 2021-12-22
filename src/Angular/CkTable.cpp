#include "Angular/CkTable.hpp"
#include "Angular/Wigner369j.hpp"
#include <cassert>
// NB: These headers MUST be in this order; otherwise fails to compile
// on macOS... what :\ At least it works now....
// Perhaps something to do with cmath vs math.h conflict???
// math.h might be included by gsl, which is in Angular_369j...
// Can only reproduce issue on mac..yay

#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace Angular {
//******************************************************************************
void CkTable::fill(const int in_max_twoj) {

  // re-factor. No longer allow 2j and k to diverge
  const int in_max_K = in_max_twoj;

  // auto max_jindex = (in_max_twoj - 1) / 2;
  auto max_jindex = Angular::jindex(in_max_twoj);
  if (max_jindex <= m_max_jindex_sofar)
    return;

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
double CkTable::get_tildeCkab_mutable(int k, int ka, int kb) {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  // parity:
  auto pi_ok = Angular::evenQ(Angular::l_k(ka) + Angular::l_k(kb) + k);
  if (!pi_ok)
    return 0;

  auto maxji = std::max(jia, jib);
  if (maxji > m_max_jindex_sofar || k > m_max_k_sofar)
    fill(std::max(k, twoj(maxji)) + 1); // XXX Hack? Why need +1 ??
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] * m_Rjab_a_b[jia][jib] :
                       m_3j_k_a_b[k][jib][jia] * m_Rjab_a_b[jib][jia];
}

double CkTable::get_tildeCkab(int k, int ka, int kb) const {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  // parity:
  auto pi_ok = Angular::evenQ(Angular::l_k(ka) + Angular::l_k(kb) + k);
  if (!pi_ok)
    return 0.0;

  if (Angular::triangle(twoj_k(ka), twoj_k(kb), 2 * k) == 0)
    return 0.0;
  assert(std::max(jia, jib) <= m_max_jindex_sofar);
  assert(k <= m_max_k_sofar);
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] * m_Rjab_a_b[jia][jib] :
                       m_3j_k_a_b[k][jib][jia] * m_Rjab_a_b[jib][jia];
}

//******************************************************************************
double CkTable::get_Ckab_mutable(int k, int ka, int kb) {
  auto s = Angular::evenQ_2(Angular::twoj_k(ka) + 1) ? 1.0 : -1.0;
  return s * get_tildeCkab_mutable(k, ka, kb);
}

double CkTable::get_Ckab(int k, int ka, int kb) const {
  auto s = Angular::evenQ_2(Angular::twoj_k(ka) + 1) ? 1.0 : -1.0;
  return s * get_tildeCkab(k, ka, kb);
}

//******************************************************************************
double CkTable::get_3jkab_mutable(int k, int ka, int kb) {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  auto maxji = std::max(jia, jib);
  if (maxji > m_max_jindex_sofar || k > m_max_k_sofar)
    fill(std::max(k, twoj(maxji)));
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] : m_3j_k_a_b[k][jib][jia];
}

double CkTable::get_3jkab(int k, int ka, int kb) const {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  assert(std::max(jia, jib) <= m_max_jindex_sofar);
  assert(k <= m_max_k_sofar);
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] : m_3j_k_a_b[k][jib][jia];
}

double CkTable::get_Lambdakab(int k, int ka, int kb) const {
  const auto pi_ok = Angular::evenQ(Angular::l_k(ka) + Angular::l_k(kb) + k);
  if (!pi_ok)
    return 0.0;
  const auto tjs = get_3jkab(k, ka, kb);
  return tjs * tjs;
}

} // namespace Angular
