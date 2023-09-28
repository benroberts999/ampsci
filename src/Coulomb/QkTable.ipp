#pragma once
#include "CoulombIntegrals.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "QkTable.hpp"
#include <algorithm>
#include <cassert>
#include <cstring> // for memcpy
#include <string_view>

namespace Coulomb {

//==============================================================================
template <Symmetry S> void CoulombTable<S>::summary() const {
  std::cout << "Summary: \n";
  int k = 0;
  auto total = 0ul;
  for (auto &qk : m_data) {
    std::cout << "k=" << k << ": " << qk.size() << " [" << qk.bucket_count()
              << "]\n";
    total += qk.size();
    ++k;
  }
  std::cout << "Total: " << total << " non-zero integrals\n";
}

template <Symmetry S> std::size_t CoulombTable<S>::count() const {
  auto total = 0ul;
  for (auto &qk : m_data) {
    total += qk.size();
  }
  return total;
}

//==============================================================================
template <Symmetry S>
void CoulombTable<S>::add(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d,
                          Real value) {
  add(k, NormalOrder(a, b, c, d), value);
}
//------------------------------------------------------------------------------
template <Symmetry S>
void CoulombTable<S>::add(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d,
                          Real value) {
  add(k, a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index(), value);
}
//------------------------------------------------------------------------------
template <Symmetry S>
void CoulombTable<S>::add(int k, nk4Index index, Real value) {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size()) {
    m_data.resize(sk + 1);
  }
  m_data.at(sk).insert({index, value});
}

//==============================================================================
template <Symmetry S>
// Updates Q in table. If not present, adds new Q
void CoulombTable<S>::update(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d,
                             Real value) {
  update(k, NormalOrder(a, b, c, d), value);
}

//------------------------------------------------------------------------------
template <Symmetry S>
// Updates Q in table. If not present, adds new Q
void CoulombTable<S>::update(int k, const DiracSpinor &a, const DiracSpinor &b,
                             const DiracSpinor &c, const DiracSpinor &d,
                             Real value) {
  update(k, a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index(), value);
}

//------------------------------------------------------------------------------
template <Symmetry S>
void CoulombTable<S>::update(int k, nk4Index index, Real value) {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size()) {
    m_data.resize(sk + 1);
  }
  m_data.at(sk).insert_or_assign(index, value);
}
//==============================================================================
template <Symmetry S>
bool CoulombTable<S>::contains(int k, nkIndex a, nkIndex b, nkIndex c,
                               nkIndex d) const {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return false;
  return m_data.at(sk).find(NormalOrder(a, b, c, d)) != m_data.at(sk).cend();
}
//------------------------------------------------------------------------------
template <Symmetry S>
bool CoulombTable<S>::contains(int k, const DiracSpinor &a,
                               const DiracSpinor &b, const DiracSpinor &c,
                               const DiracSpinor &d) const {
  return contains(k, a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
}
//------------------------------------------------------------------------------
template <Symmetry S>
bool CoulombTable<S>::contains(int k, nk4Index index) const {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return false;
  return m_data.at(sk).find(index) != m_data.at(sk).cend();
}

//==============================================================================
//==============================================================================

template <Symmetry S> double *CoulombTable<S>::get(int k, nk4Index index) {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return nullptr;
  // check valid k? Probably faster to lookup in table
  const auto map_it = m_data.at(sk).find(index);
  if (map_it == m_data.at(sk).cend())
    return nullptr;
  return &(map_it->second);
}

//==============================================================================
template <Symmetry S>
double CoulombTable<S>::Q(int k, nkIndex a, nkIndex b, nkIndex c,
                          nkIndex d) const {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return 0.0;
  // check valid k? Probably faster to lookup in table
  const auto map_it = m_data.at(sk).find(NormalOrder(a, b, c, d));
  if (map_it == m_data.at(sk).cend())
    return 0.0;
  return map_it->second;
}
//------------------------------------------------------------------------------
template <Symmetry S>
double CoulombTable<S>::Q(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d) const {
  return Q(k, a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
}
//------------------------------------------------------------------------------
template <Symmetry S> double CoulombTable<S>::Q(int k, nk4Index index) const {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return 0.0;
  // check valid k? Probably faster to lookup in table
  const auto map_it = m_data.at(sk).find(index);
  if (map_it == m_data.at(sk).cend())
    return 0.0;
  return map_it->second;
}

//==============================================================================
template <Symmetry S>
double CoulombTable<S>::R(int k, nkIndex a, nkIndex b, nkIndex c,
                          nkIndex d) const {
  const auto tQk = Q(k, a, b, c, d);
  if (tQk == 0.0)
    return 0.0;
  const auto s = Angular::neg1pow(k);

  const auto [na, ka] = Angular::index_to_nk(a);
  const auto [nb, kb] = Angular::index_to_nk(b);
  const auto [nc, kc] = Angular::index_to_nk(c);
  const auto [nd, kd] = Angular::index_to_nk(d);

  const auto tCkac = Angular::tildeCk_kk(k, ka, kc);
  const auto tCkbd = Angular::tildeCk_kk(k, kb, kd);
  return tQk / (s * tCkac * tCkbd);
}
//------------------------------------------------------------------------------
template <Symmetry S>
double CoulombTable<S>::R(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d) const {
  return R(k, a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
}

//==============================================================================
template <Symmetry S>
double CoulombTable<S>::P(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d,
                          const Angular::SixJTable *const sj) const {
  double Pk_abcd{0.0};

  const auto ka = Angular::nkindex_to_kappa(a);
  const auto kb = Angular::nkindex_to_kappa(b);
  const auto kc = Angular::nkindex_to_kappa(c);
  const auto kd = Angular::nkindex_to_kappa(d);

  const auto tja = Angular::twoj_k(ka);
  const auto tjb = Angular::twoj_k(kb);
  const auto tjc = Angular::twoj_k(kc);
  const auto tjd = Angular::twoj_k(kd);

  // 6j(s) Triads: {a,c,k}, {k,b,d}, {c,b,l}, {d,a,l}
  if (Angular::triangle(tja, tjc, 2 * k) == 0 ||
      Angular::triangle(2 * k, tjb, tjd) == 0)
    return 0.0;

  const auto [lmin, lmax] = k_minmax_Q(ka, kb, kd, kc); // exchange
  for (int l = lmin; l <= lmax; l += 2) {
    const auto ql = this->Q(l, a, b, d, c); // exchange
    if (ql == 0.0)
      continue;
    const auto sixj = sj ? sj->get_2(tja, tjc, 2 * k, tjb, tjd, 2 * l) :
                           Angular::sixj_2(tja, tjc, 2 * k, tjb, tjd, 2 * l);
    Pk_abcd += sixj * ql;
  }
  Pk_abcd *= double(2 * k + 1);
  return Pk_abcd;
}

//------------------------------------------------------------------------------

template <Symmetry S>
double CoulombTable<S>::P(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d,
                          const Angular::SixJTable *const sj) const {

  return P(k, a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index(), sj);
}

//==============================================================================
template <Symmetry S>
double CoulombTable<S>::P2(int k, const DiracSpinor &a, const DiracSpinor &b,
                           const DiracSpinor &c, const DiracSpinor &d,
                           const Angular::SixJTable &sj,
                           const std::vector<double> &fk) const {
  double Pk_abcd{0.0};

  // 6j(s) Triads: {a,c,k}, {k,b,d}, {c,b,l}, {d,a,l}
  if (Coulomb::triangle(a, c, k) == 0 || Coulomb::triangle(k, b, d) == 0)
    return 0.0;

  const auto [lmin, lmax] = k_minmax_Q(a, b, d, c); // exchange
  for (int l = lmin; l <= lmax; l += 2) {
    const auto ql = this->Q(l, a, b, d, c); // exchange
    if (ql == 0.0)
      continue;
    const auto sixj = sj.get(a, c, k, b, d, l);

    // include effective Coulomb screening:
    const auto f_scr_l = (l < (int)fk.size()) ? fk[std::size_t(l)] : 1.0;

    Pk_abcd += f_scr_l * sixj * ql;
  }
  Pk_abcd *= double(2 * k + 1);
  return Pk_abcd;
}

//==============================================================================
template <Symmetry S>
double CoulombTable<S>::W(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d,
                          const Angular::SixJTable *const sj) const {
  return W(k, a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index(), sj);
}

template <Symmetry S>
double CoulombTable<S>::W(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d,
                          const Angular::SixJTable *const sj) const {
  return Q(k, a, b, c, d) + P(k, a, b, c, d, sj);
}

//==============================================================================
template <Symmetry S>
double CoulombTable<S>::g(const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d, int tma,
                          int tmb, int tmc, int tmd) const {

  if (tmc - tma != tmb - tmd)
    return 0.0;
  const int twoq = tmc - tma;
  const auto s = Angular::neg1pow_2(tma - tmb + twoq);

  //
  double g = 0.0;
  const auto [k0, ki] = k_minmax_Q(a, b, c, d);
  for (int k = k0; k <= ki; k += 2) {
    const auto tjs1 =
        Angular::threej_2(a.twoj(), 2 * k, c.twoj(), -tma, -twoq, tmc);
    const auto tjs2 =
        Angular::threej_2(b.twoj(), 2 * k, d.twoj(), -tmb, twoq, tmd);
    if (tjs1 == 0.0 || tjs2 == 0.0)
      continue;
    g += Angular::neg1pow(k) * s * tjs1 * tjs2 * Q(k, a, b, c, d);
  }
  return g;
}

//==============================================================================
//==============================================================================

//==============================================================================
template <>
inline nk4Index
CoulombTable<Symmetry::Qk>::NormalOrder_impl(nkIndex a, nkIndex b, nkIndex c,
                                             nkIndex d) {

  // abcd -> ijkl, with i = min(a,b,c,d)

  // abcd = adcb
  // badc = bcda
  // cbad = cdab
  // dabc = dcba
  // nb: there must be a more efficient way of doing this!

  // const auto tmp1 = FormIndex(a, b, c, d);
  // const auto tmp2 = FormIndex(a, d, c, b);
  // const auto tmp3 = FormIndex(b, a, d, c);
  // const auto tmp4 = FormIndex(b, c, d, a);
  // const auto tmp5 = FormIndex(c, b, a, d);
  // const auto tmp6 = FormIndex(c, d, a, b);
  // const auto tmp7 = FormIndex(d, a, b, c);
  // const auto tmp8 = FormIndex(d, c, b, a);
  // return std::min({tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8});

  nk4Index out = FormIndex(a, b, c, d);
  const auto min = std::min({a, b, c, d});
  if (a == min) {
    const auto tmp = FormIndex(a, d, c, b);
    if (tmp < out)
      out = tmp;
  }
  if (b == min) {
    const auto tmp = std::min(FormIndex(b, a, d, c), FormIndex(b, c, d, a));
    if (tmp < out)
      out = tmp;
  }
  if (c == min) {
    const auto tmp = std::min(FormIndex(c, b, a, d), FormIndex(c, d, a, b));
    if (tmp < out)
      out = tmp;
  }
  if (d == min) {
    const auto tmp = std::min(FormIndex(d, a, b, c), FormIndex(d, c, b, a));
    if (tmp < out)
      out = tmp;
  }
  return out;
}
//------------------------------------------------------------------------------
template <>
inline nk4Index
CoulombTable<Symmetry::Wk>::NormalOrder_impl(nkIndex a, nkIndex b, nkIndex c,
                                             nkIndex d) {
  const auto tmp1 = FormIndex(a, b, c, d);
  const auto tmp2 = FormIndex(b, a, d, c);
  const auto tmp3 = FormIndex(c, d, a, b);
  const auto tmp4 = FormIndex(d, c, b, a);
  return std::min({tmp1, tmp2, tmp3, tmp4});
}

//------------------------------------------------------------------------------
template <>
inline nk4Index
CoulombTable<Symmetry::Lk>::NormalOrder_impl(nkIndex a, nkIndex b, nkIndex c,
                                             nkIndex d) {
  const auto tmp1 = FormIndex(a, b, c, d);
  const auto tmp2 = FormIndex(b, a, d, c);
  return std::min({tmp1, tmp2});
}

//------------------------------------------------------------------------------
template <>
inline nk4Index
CoulombTable<Symmetry::none>::NormalOrder_impl(nkIndex a, nkIndex b, nkIndex c,
                                               nkIndex d) {
  return FormIndex(a, b, c, d);
}

//==============================================================================
template <Symmetry S>
nk4Index CoulombTable<S>::NormalOrder(nkIndex a, nkIndex b, nkIndex c,
                                      nkIndex d) const {
  return NormalOrder_impl(a, b, c, d);
}
//------------------------------------------------------------------------------
template <Symmetry S>
nk4Index
CoulombTable<S>::NormalOrder(const DiracSpinor &a, const DiracSpinor &b,
                             const DiracSpinor &c, const DiracSpinor &d) const {
  return NormalOrder(a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
}

//==============================================================================
template <Symmetry S>
nk4Index CoulombTable<S>::CurrentOrder(nkIndex a, nkIndex b, nkIndex c,
                                       nkIndex d) const {
  return FormIndex(a, b, c, d);
}
//------------------------------------------------------------------------------
template <Symmetry S>
nk4Index CoulombTable<S>::CurrentOrder(const DiracSpinor &a,
                                       const DiracSpinor &b,
                                       const DiracSpinor &c,
                                       const DiracSpinor &d) const {
  return CurrentOrder(a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
}

//==============================================================================
template <Symmetry S>
nk4Index CoulombTable<S>::FormIndex(nkIndex a, nkIndex b, nkIndex c,
                                    nkIndex d) {
  // create a nk4Index from four small nkIndex, by laying out in memory

  // static_assert(sizeof(nk4Index) == 4 * sizeof(nkIndex));
  // nk4Index big_index;
  // const auto p_index = reinterpret_cast<nkIndex *>(&big_index);
  // std::memcpy(p_index, &d, sizeof(nkIndex));
  // std::memcpy(p_index + 1, &c, sizeof(nkIndex));
  // std::memcpy(p_index + 2, &b, sizeof(nkIndex));
  // std::memcpy(p_index + 3, &a, sizeof(nkIndex));
  // return big_index;

  // this seems slightly faster...though variance large
  static_assert(sizeof(nk4Index) == 4 * sizeof(nkIndex));
  static_assert(sizeof(nkIndex) * 8 == 16);
  return ((nk4Index)a << 48) + ((nk4Index)b << 32) + ((nk4Index)c << 16) +
         (nk4Index)d;
}

//==============================================================================
template <Symmetry S>
std::array<nkIndex, 4>
CoulombTable<S>::UnFormIndex(const nk4Index &index) const {
  static_assert(sizeof(nk4Index) == sizeof(std::array<nkIndex, 4>));
  std::array<nkIndex, 4> set;
  std::memcpy(&set, &index, sizeof(nk4Index));
  return set;
}

//==============================================================================
//==============================================================================
template <Symmetry S>
void CoulombTable<S>::fill(const std::vector<DiracSpinor> &basis,
                           const YkTable &yk, int k_cut, bool print) {
  static_assert(S == Symmetry::Qk);
  IO::ChronoTimer t("fill");

  /*
  In order to make best use of CPU and Memory, we "fill" the QK table in a
  strange 4-step manner.
  1. Count the number of non-zero Q's for each k
  2. Use info from aboe to 'reserve()' space in the map
  3. Fill all the non-zero Q entries in the map with 0. Note that this adds new
  elements into the map, so cannot be done in parallel. Since we store a
  different map for each k, we may parallelise over k (though this is not super
  efficient)
  4. Now that we have a fully sized map, we can update each of its values in
  parallel
  Note that we make use of the symmetry to ensure we do not access any element
  from more than 1 thread (also, allowing us to not do any more calculations
  than necisary)
  */

  const auto tmp_max_k =
      std::size_t(std::max(DiracSpinor::max_tj(basis), 1) - 1);

  const auto max_k =
      (k_cut <= 0) ? tmp_max_k : std::min(tmp_max_k, std::size_t(k_cut));

  if (m_data.size() < max_k + 1)
    m_data.resize(max_k + 1);

  // 1) Count non-zero Q integrals (each k). Use this to 'reserve' map space
  t.start();
  std::vector<std::size_t> count_non_zero_k(max_k + 1);
#pragma omp parallel for
  for (auto k = 0ul; k <= max_k; ++k) {
    const auto ik = static_cast<int>(k);
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        for (const auto &c : basis) {
          if (!Angular::Ck_kk_SR(ik, a.kappa(), c.kappa()))
            continue;
          for (const auto &d : basis) {
            // due to symmetry, only calculate each 'unique' integral once
            if (NormalOrder(a, b, c, d) == CurrentOrder(a, b, c, d)) {
              if (Angular::Ck_kk_SR(ik, b.kappa(), d.kappa())) {
                ++count_non_zero_k[k];
              }
            }
          }
        }
      }
    }
  }
  if (print)
    std::cout << "Count non-zero: " << t.lap_reading_str() << std::endl;

  // 2) Reserve space in each sub-map
  t.start();
#pragma omp parallel for
  for (auto ik = 0ul; ik <= max_k; ++ik) {
    m_data[ik].reserve(count_non_zero_k[ik]);
  }
  if (print)
    std::cout << "Reserve: " << t.lap_reading_str() << std::endl;

  // 3) Create space in map (set each element to zero).
  t.start();
#pragma omp parallel for
  for (auto k = 0ul; k <= max_k; ++k) {
    for (const auto &a : basis) {
      for (const auto &c : basis) {
        if (!Angular::Ck_kk_SR(int(k), a.kappa(), c.kappa()))
          continue;
        for (const auto &b : basis) {
          for (const auto &d : basis) {
            if (!Angular::Ck_kk_SR(int(k), b.kappa(), d.kappa()))
              continue;
            const auto normal_index = NormalOrder(a, b, c, d);
            if (normal_index == CurrentOrder(a, b, c, d)) {
              add(int(k), normal_index, 0.0);
            }
          }
        }
      }
    }
  }
  if (print)
    std::cout << "Fill w/ zeros: " << t.lap_reading_str() << std::endl;

  // 4) Fill the pre-constructed map with values, in parallel. Since we are
  // not adding any new elements to map, and since we are guarenteed to only
  // access each map element once, we can do this part in parallel. nb: This
  // //isation is not very efficient, though in theory it can be 100%
  if (print)
    std::cout << "Fill w/ values: " << std::flush;
  t.start();
#pragma omp parallel for collapse(2)
  for (auto ia = 0ul; ia < basis.size(); ++ia) {
    for (auto ib = 0ul; ib < basis.size(); ++ib) {
      const auto &a = basis[ia];
      const auto &b = basis[ib];
      for (const auto &c : basis) {
        for (const auto &d : basis) {
          const auto normal_index = NormalOrder(a, b, c, d);
          if (normal_index == CurrentOrder(a, b, c, d)) {
            const auto [kmin, kmax] = k_minmax_Q(a, b, c, d);
            for (int k = kmin; k <= kmax && k <= int(max_k); k += 2) {
              double *ptr = get(k, normal_index);
              assert(ptr != nullptr);
              if (*ptr == 0.0) {
                // only calculate if not already in table
                // This saves some time, but surprisingly little!
                *ptr = yk.Q(k, a, b, c, d);
              }
              // update(k, a, b, c, d, yk.Q(k, a, b, c, d));
            }
          }
        }
      }
    }
  }

  if (print)
    std::cout << t.lap_reading_str() << std::endl;

  if (print)
    summary();
}

//==============================================================================
template <Symmetry S>
void CoulombTable<S>::fill(const std::vector<DiracSpinor> &basis,
                           const CoulombFunction &Fk,
                           const SelectionRules &Fk_SR, int k_cut, bool print) {
  IO::ChronoTimer t("fill");

  /*
  In order to make best use of CPU and Memory, we "fill" the QK table in a
  strange 4-step manner.
  1. Count the number of non-zero Q's for each k
  2. Use info from aboe to 'reserve()' space in the map
  3. Fill all the non-zero Q entries in the map with 0. Note that this adds new
  elements into the map, so cannot be done in parallel. Since we store a
  different map for each k, we may parallelise over k (though this is not super
  efficient)
  4. Now that we have a fully sized map, we can update each of its values in
  parallel
  Note that we make use of the symmetry to ensure we do not access any element
  from more than 1 thread (also, allowing us to not do any more calculations
  than necisary)
  */

  const auto tmp_max_k = std::size_t(DiracSpinor::max_tj(basis));
  // May have different parity rule, so don't -1 here

  const auto max_k =
      (k_cut <= 0) ? tmp_max_k : std::min(tmp_max_k, std::size_t(k_cut));

  if (m_data.size() < max_k + 1)
    m_data.resize(max_k + 1);

  // 1) Count non-zero Q integrals (each k). Use this to 'reserve' map space
  t.start();
  std::vector<std::size_t> count_non_zero_k(max_k + 1);
#pragma omp parallel for
  for (auto k = 0ul; k <= max_k; ++k) {
    const auto ik = static_cast<int>(k);
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        for (const auto &c : basis) {
          for (const auto &d : basis) {
            // due to symmetry, only calculate each 'unique' integral once
            if (NormalOrder(a, b, c, d) == CurrentOrder(a, b, c, d)) {
              if (Fk_SR(ik, a, b, c, d)) {
                ++count_non_zero_k[k];
              }
            }
          }
        }
      }
    }
  }
  if (print)
    std::cout << "Count non-zero: " << t.lap_reading_str() << std::endl;

  // 2) Reserve space in each sub-map
  t.start();
#pragma omp parallel for
  for (auto ik = 0ul; ik <= max_k; ++ik) {
    m_data[ik].reserve(count_non_zero_k[ik]);
  }
  if (print)
    std::cout << "Reserve: " << t.lap_reading_str() << std::endl;

  // 3) Create space in map (set each element to zero).
  t.start();
#pragma omp parallel for
  for (auto k = 0ul; k <= max_k; ++k) {
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        for (const auto &c : basis) {
          for (const auto &d : basis) {
            if (Fk_SR(int(k), a, b, c, d)) {
              const auto normal_index = NormalOrder(a, b, c, d);
              if (normal_index == CurrentOrder(a, b, c, d)) {
                add(int(k), normal_index, 0.0);
              }
            }
          }
        }
      }
    }
  }
  if (print)
    std::cout << "Fill w/ zeros: " << t.lap_reading_str() << std::endl;

  // 4) Fill the pre-constructed map with values, in parallel. Since we are
  // not adding any new elements to map, and since we are guarenteed to only
  // access each map element once, we can do this part in parallel. nb: This
  // //isation is not very efficient, though in theory it can be 100%
  if (print)
    std::cout << "Fill w/ values: " << std::flush;
  t.start();

#pragma omp parallel for collapse(2)
  for (const auto &a : basis) {
    // const auto &a = basis[ia];
    // t.start();
    // For tests only:
    // std::cout << ia << "/" << basis.size() << "     \r" << std::flush;
    // Faster to parelise here??? or random?
    for (const auto &b : basis) {
      for (const auto &c : basis) {
        for (const auto &d : basis) {
          const auto normal_index = NormalOrder(a, b, c, d);
          if (normal_index == CurrentOrder(a, b, c, d)) {
            for (int k = 0; k <= int(max_k); ++k) {
              if (Fk_SR(k, a, b, c, d)) {
                double *ptr = get(k, normal_index);
                assert(ptr != nullptr);
                if (*ptr == 0.0) {
                  // only calculate if not already in table
                  *ptr = Fk(k, a, b, c, d);
                }
                // update(k, a, b, c, d, Fk(k, a, b, c, d));
              }
            }
          }
        }
      }
    }
  }
  if (print)
    std::cout << t.lap_reading_str() << std::endl;

  if (print)
    summary();
}

//==============================================================================
template <Symmetry S>
void CoulombTable<S>::write(const std::string &fname) const {
  if (fname == "false")
    return;
  std::cout << "Writing " << count() << " integrals to file: " << fname << ".."
            << std::flush;
  std::fstream f;
  const auto rw = IO::FRW::write;
  IO::FRW::open_binary(f, fname, rw);

  auto size = m_data.size();
  rw_binary(f, rw, size);
  for (const auto &Q_k : m_data) {
    auto size_k = Q_k.size();
    rw_binary(f, rw, size_k);
    for (auto [key, value] : Q_k) {
      auto key_copy = key; // have no pass non-const reference!
      rw_binary(f, rw, key_copy, value);
    }
  }
  std::cout << "\n" << std::flush;
}

//==============================================================================
template <Symmetry S> bool CoulombTable<S>::read(const std::string &fname) {
  if (fname == "false")
    return false;
  std::fstream f;
  const auto rw = IO::FRW::read;
  IO::FRW::open_binary(f, fname, rw);

  if (!f.good())
    return false;

  std::size_t size{0};
  rw_binary(f, rw, size);
  m_data.resize(size);
  if (m_data.size() == 0)
    return false;
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    auto &Q_k = m_data[i];
    std::size_t size_k{0};
    rw_binary(f, rw, size_k);
    Q_k.reserve(size_k);
    for (std::size_t ik = 0; ik < size_k; ++ik) {
      nk4Index key;
      Real value;
      rw_binary(f, rw, key, value);
      Q_k[key] = value;
    }
  }
  std::cout << "Read " << count() << " integrals from file: " << fname << "\n"
            << std::flush;
  return true;
}

} // namespace Coulomb
