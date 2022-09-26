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
  std::cout << "total: " << total << " non-zero Qk\n";
}

template <Symmetry S> int CoulombTable<S>::count() const {
  auto total = 0ul;
  for (auto &qk : m_data) {
    total += qk.size();
  }
  return static_cast<int>(total);
}

//==============================================================================
template <Symmetry S>
void CoulombTable<S>::add(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d,
                          Real value) {
  add(k, NormalOrder(a, b, c, d), value);
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
void CoulombTable<S>::update(int k, const DiracSpinor &a, const DiracSpinor &b,
                             const DiracSpinor &c, const DiracSpinor &d,
                             Real value) {
  update(k, NormalOrder(a, b, c, d), value);
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
bool CoulombTable<S>::contains(int k, const DiracSpinor &a,
                               const DiracSpinor &b, const DiracSpinor &c,
                               const DiracSpinor &d) const {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return false;
  return m_data.at(sk).find(NormalOrder(a, b, c, d)) != m_data.at(sk).cend();
}

template <Symmetry S>
bool CoulombTable<S>::contains(int k, nk4Index index) const {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return false;
  return m_data.at(sk).find(index) != m_data.at(sk).cend();
}

//==============================================================================
//==============================================================================

//==============================================================================
template <Symmetry S>
double CoulombTable<S>::Q(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d) const {
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
double CoulombTable<S>::R(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d) const {
  const auto tQk = Q(k, a, b, c, d);
  if (tQk == 0.0)
    return 0.0;
  const auto s = Angular::neg1pow(k);
  const auto tCkac = Angular::tildeCk_kk(k, a.kappa(), c.kappa());
  const auto tCkbd = Angular::tildeCk_kk(k, b.kappa(), d.kappa());
  return tQk / (s * tCkac * tCkbd);
}

//==============================================================================
template <Symmetry S>
double CoulombTable<S>::P(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d,
                          const Angular::SixJTable *const sj) const {
  double Pk_abcd{0.0};

  // 6j(s) Triads: {a,c,k}, {k,b,d}, {c,b,l}, {d,a,l}
  if (Coulomb::triangle(a, c, k) == 0 || Coulomb::triangle(k, b, d) == 0)
    return 0.0;

  const auto [lmin, lmax] = k_minmax_Q(a, b, d, c); // exchange
  for (int l = lmin; l <= lmax; l += 2) {
    const auto ql = this->Q(l, a, b, d, c); // exchange
    if (ql == 0.0)
      continue;
    const auto sixj =
        sj ? sj->get(a, c, k, b, d, l) : Coulomb::sixj(a, c, k, b, d, l);
    Pk_abcd += sixj * ql;
  }
  Pk_abcd *= double(2 * k + 1);
  return Pk_abcd;
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
  return Q(k, a, b, c, d) + P(k, a, b, c, d, sj);
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
nk4Index
CoulombTable<S>::NormalOrder(const DiracSpinor &a, const DiracSpinor &b,
                             const DiracSpinor &c, const DiracSpinor &d) const {
  return NormalOrder_impl(a.nk_index(), b.nk_index(), c.nk_index(),
                          d.nk_index());
  // if constexpr (S == Coulomb::Symmetry::Qk)
  //   return NormalOrder_Q(a.nk_index(), b.nk_index(), c.nk_index(),
  //                        d.nk_index());
  // if constexpr (S == Coulomb::Symmetry::Wk)
  //   return NormalOrder_W(a.nk_index(), b.nk_index(), c.nk_index(),
  //                        d.nk_index());
  // if constexpr (S == Coulomb::Symmetry::Lk)
  //   return NormalOrder_L(a.nk_index(), b.nk_index(), c.nk_index(),
  //                        d.nk_index());
  // if constexpr (S == Coulomb::Symmetry::none)
  //   return NormalOrder_none(a.nk_index(), b.nk_index(), c.nk_index(),
  //                           d.nk_index());
}

//==============================================================================
template <Symmetry S>
nk4Index CoulombTable<S>::CurrentOrder(const DiracSpinor &a,
                                       const DiracSpinor &b,
                                       const DiracSpinor &c,
                                       const DiracSpinor &d) const {
  return FormIndex(a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
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
  return (nk4Index)d + ((nk4Index)c << 16) + ((nk4Index)b << 32) +
         ((nk4Index)a << 48);
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
                           const YkTable &yk, int k_cut) {
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

  const auto tmp_max_k = std::size_t(DiracSpinor::max_tj(basis) - 1);

  const auto max_k =
      (k_cut <= 0) ? tmp_max_k : std::min(tmp_max_k, std::size_t(k_cut));

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
  std::cout << "Count non-zero: " << t.lap_reading_str() << std::endl;

  // 2) Reserve space in each sub-map
  t.start();
#pragma omp parallel for
  for (auto ik = 0ul; ik <= max_k; ++ik) {
    m_data[ik].reserve(count_non_zero_k[ik]);
  }
  std::cout << "Reserve: " << t.lap_reading_str() << std::endl;

  // 3) Create space in map (set each element to zero).
  t.start();
#pragma omp parallel for
  for (auto k = 0ul; k <= max_k; ++k) {
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        for (const auto &c : basis) {
          if (!Angular::Ck_kk_SR(int(k), a.kappa(), c.kappa()))
            continue;
          for (const auto &d : basis) {
            if (NormalOrder(a, b, c, d) == CurrentOrder(a, b, c, d)) {
              if (Angular::Ck_kk_SR(int(k), b.kappa(), d.kappa())) {
                add(int(k), a, b, c, d, 0.0);
              }
            }
          }
        }
      }
    }
  }
  std::cout << "Fill w/ zeros: " << t.lap_reading_str() << std::endl;

  // 4) Fill the pre-constructed map with values, in parallel. Since we are
  // not adding any new elements to map, and since we are guarenteed to only
  // access each map element once, we can do this part in parallel. nb: This
  // //isation is not very efficient, though in theory it can be 100%
  t.start();
#pragma omp parallel for
  for (auto ia = 0ul; ia < basis.size(); ++ia) {
    const auto &a = basis[ia];
    for (const auto &b : basis) {
      for (const auto &c : basis) {
        for (const auto &d : basis) {
          if (NormalOrder(a, b, c, d) == CurrentOrder(a, b, c, d)) {
            auto [kmin, kmax] = k_minmax_Q(a, b, c, d);
            kmax = std::clamp(kmax, 0, int(max_k));
            for (int k = kmin; k <= kmax; k += 2) {
              update(k, a, b, c, d, yk.Q(k, a, b, c, d));
            }
          }
        }
      }
    }
  }

  std::cout << "Fill w/ values: " << t.lap_reading_str() << std::endl;

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
