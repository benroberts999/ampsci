#include "Physics/Angular.hpp"

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
        m_3j_k_a_b[k][jia].push_back(Wigner::special_threej_2(tja, tjb, 2 * k));
      }
    }
  }

  m_max_jindex_sofar = max_jindex;
  m_max_k_sofar = in_max_K;
}

//******************************************************************************
double Ck_ab::get_tildeCkab(int k, int ka, int kb) {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  // parity:
  auto pi_ok = Wigner::evenQ(Wigner::l_k(ka) + Wigner::l_k(kb) + k);
  if (!pi_ok)
    return 0;

  auto maxji = std::max(jia, jib);
  if (maxji > m_max_jindex_sofar || k > m_max_k_sofar)
    fill_maxK_twojmax(k, 2 * maxji + 1);
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] * m_Rjab_a_b[jia][jib]
                     : m_3j_k_a_b[k][jib][jia] * m_Rjab_a_b[jib][jia];
}

double Ck_ab::get_tildeCkab_const(int k, int ka, int kb) const {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  // parity:
  auto pi_ok = Wigner::evenQ(Wigner::l_k(ka) + Wigner::l_k(kb) + k);
  if (!pi_ok)
    return 0;
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] * m_Rjab_a_b[jia][jib]
                     : m_3j_k_a_b[k][jib][jia] * m_Rjab_a_b[jib][jia];
}

//******************************************************************************
double Ck_ab::get_Ckab(int k, int ka, int kb) {
  auto s = Wigner::evenQ_2(Wigner::twoj_k(ka) + 1) ? 1.0 : -1.0;
  return s * get_tildeCkab(k, ka, kb);
}

double Ck_ab::get_Ckab_const(int k, int ka, int kb) const {
  auto s = Wigner::evenQ_2(Wigner::twoj_k(ka) + 1) ? 1.0 : -1.0;
  return s * get_tildeCkab_const(k, ka, kb);
}

//******************************************************************************
double Ck_ab::get_3jkab(int k, int ka, int kb) {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  auto maxji = std::max(jia, jib);
  if (maxji > m_max_jindex_sofar || k > m_max_k_sofar)
    fill_maxK_twojmax(k, 2 * maxji + 1);
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] : m_3j_k_a_b[k][jib][jia];
}

double Ck_ab::get_3jkab_const(int k, int ka, int kb) const {
  auto jia = jindex_kappa(ka);
  auto jib = jindex_kappa(kb);
  return (jia > jib) ? m_3j_k_a_b[k][jia][jib] : m_3j_k_a_b[k][jib][jia];
}

} // namespace Angular
