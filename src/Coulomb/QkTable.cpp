#include "QkTable.hpp"
#include "CoulombIntegrals.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include <algorithm>
#include <cassert>
#include <cstring> // for memcpy
#include <string_view>

namespace Coulomb {

//******************************************************************************
void CoulombTable::count() const {
  std::cout << "Count: \n";
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

//******************************************************************************
void CoulombTable::add(int k, const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d, Real value) {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size()) {
    m_data.resize(sk + 1);
  }
  m_data.at(sk).insert({NormalOrder(a, b, c, d), value});
}
//----------------
void CoulombTable::add(int k, BigIndex index, Real value) {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size()) {
    m_data.resize(sk + 1);
  }
  m_data.at(sk).insert({index, value});
}

//******************************************************************************
// Updates Q in table. If not present, adds new Q
void CoulombTable::update(int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d,
                          Real value) {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size()) {
    m_data.resize(sk + 1);
  }
  m_data.at(sk).insert_or_assign(NormalOrder(a, b, c, d), value);
}

//******************************************************************************
bool CoulombTable::contains(int k, const DiracSpinor &a, const DiracSpinor &b,
                            const DiracSpinor &c, const DiracSpinor &d) const {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return false;
  return m_data.at(sk).find(NormalOrder(a, b, c, d)) != m_data.at(sk).cend();
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
double CoulombTable::Q(int k, const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d) const {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return 0.0;
  // check valid k? Probably faster to lookup in table
  const auto map_it = m_data.at(sk).find(NormalOrder(a, b, c, d));
  if (map_it == m_data.at(sk).end())
    return 0.0;
  return map_it->second;
}

//******************************************************************************
double CoulombTable::R(int k, const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d) const {
  const auto tQk = Q(k, a, b, c, d);
  if (tQk == 0.0)
    return 0.0;
  const auto s = Angular::neg1pow(k);
  const auto tCkac = Angular::tildeCk_kk(k, a.k, c.k);
  const auto tCkbd = Angular::tildeCk_kk(k, b.k, d.k);
  return tQk / (s * tCkac * tCkbd);
}

//******************************************************************************
double CoulombTable::P(int k, const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d) const {
  double Pk_abcd{0.0};
  const auto [lmin, lmax] = k_minmax_Q(a, b, d, c); // exchange
  for (int l = lmin; l <= lmax; l += 2) {
    const auto ql = Q(l, a, b, d, c); // exchange
    if (ql == 0.0)
      continue;
    const auto sixj = Coulomb::sixj(a, c, k, b, d, l);
    Pk_abcd += sixj * ql;
  }
  Pk_abcd *= double(2 * k + 1);
  return Pk_abcd;
}

//******************************************************************************
double CoulombTable::W(int k, const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d) const {
  return Q(k, a, b, c, d) + P(k, a, b, c, d);
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
CoulombTable::BigIndex CoulombTable::NormalOrder(const DiracSpinor &a,
                                                 const DiracSpinor &b,
                                                 const DiracSpinor &c,
                                                 const DiracSpinor &d) const {
  // xxx conversion!
  return NormalOrder_impl(Index(a.nk_index()), Index(b.nk_index()),
                          Index(c.nk_index()), Index(d.nk_index()));
}

//******************************************************************************
CoulombTable::BigIndex CoulombTable::CurrentOrder(const DiracSpinor &a,
                                                  const DiracSpinor &b,
                                                  const DiracSpinor &c,
                                                  const DiracSpinor &d) const {
  // xxx conversion!
  return FormIndex({Index(a.nk_index()), Index(b.nk_index()),
                    Index(c.nk_index()), Index(d.nk_index())});
}

//******************************************************************************
CoulombTable::BigIndex QkTable::NormalOrder_impl(Index a, Index b, Index c,
                                                 Index d) const {
  // put smallest first
  const auto min = std::min({a, b, c, d});
  if (min == a) {
    // options are abcd, and adcb
    return (b < d) ? FormIndex({a, b, c, d}) : FormIndex({a, d, c, b});
  } else if (min == b) {
    // options are badc, and bcda
    return (a < c) ? FormIndex({b, a, d, c}) : FormIndex({b, c, d, a});
  } else if (min == c) {
    // options are cbad, and cdab
    return (b < d) ? FormIndex({c, b, a, d}) : FormIndex({c, d, a, b});
  } else if (min == d) {
    // options are dabc, and dcba
    return (a < c) ? FormIndex({d, a, b, c}) : FormIndex({d, c, b, a});
  }
  assert(false);
}
//------------------------------------------------------------------------------
CoulombTable::BigIndex WkTable::NormalOrder_impl(Index a, Index b, Index c,
                                                 Index d) const {
  // put smallest first
  const auto min = std::min({a, b, c, d});
  if (min == a) {
    return FormIndex({a, b, c, d});
  } else if (min == b) {
    return FormIndex({b, a, d, c});
  } else if (min == c) {
    return FormIndex({c, d, a, b});
  } else if (min == d) {
    return FormIndex({d, c, b, a});
  }
  assert(false);
}

//------------------------------------------------------------------------------
CoulombTable::BigIndex NkTable::NormalOrder_impl(Index a, Index b, Index c,
                                                 Index d) const {
  return FormIndex({a, b, c, d});
}

//******************************************************************************
CoulombTable::BigIndex CoulombTable::FormIndex(const IndexSet &set) const { //
  static_assert(sizeof(BigIndex) == 4 * sizeof(Index));
  static_assert(sizeof(BigIndex) == sizeof(IndexSet));
  // nb: 2x slower than necisary...forms array, _then_ copies memory..
  BigIndex index;
  // void* memcpy( void* dest, const void* src, std::size_t count )
  std::memcpy(&index, &set, sizeof(BigIndex));
  return index;
}

//******************************************************************************
CoulombTable::IndexSet CoulombTable::UnFormIndex(const BigIndex &index) const {
  IndexSet set;
  // void* memcpy( void* dest, const void* src, std::size_t count )
  std::memcpy(&set, &index, sizeof(BigIndex));
  return set;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
void CoulombTable::fill(const std::vector<DiracSpinor> &basis,
                        const YkTable &yk) {
  IO::ChronoTimer t("fill");
  // assert(yk.a_is_b());
  // const auto &basis = yk.get_a();

  // XXX Note: This uses 2x memory.. issue?
  // but, much faster than not!

  // Allow fill in parallel, by first storing in a vector, then adding vector
  // to map in series
  using TMP = std::pair<BigIndex, Real>;
  std::vector<std::vector<std::vector<TMP>>> maps_k_a;

  const auto max_k = std::size_t(DiracSpinor::max_tj(basis));
  maps_k_a.resize(max_k + 1);
  for (auto &mk : maps_k_a) {
    mk.resize(basis.size());
  }

  // Cannot insert into map in thread-safe manner.
  // So, fill a std::vector safely, then insert those elements into map
  // In order to avoid calculating equivilant Qk's twice (due to symmetry), only
  // calculate when already in NormalOrder

#pragma omp parallel for
  for (auto i = 0ul; i < basis.size(); ++i) {
    const auto &a = basis[i];
    for (const auto &b : basis) {
      for (const auto &c : basis) {
        for (const auto &d : basis) {

          // enfore symmetry here, to avoid calculating anything twice
          const auto ix = NormalOrder(a, b, c, d);
          if (ix != CurrentOrder(a, b, c, d))
            continue;

          const auto [kmin, kmax] = k_minmax_Q(a, b, c, d);
          for (int k = kmin; k <= kmax; k += 2) {
            const auto yk_bd = yk.get(k, b, d);
            if (yk_bd == nullptr)
              continue;

            const auto qk = yk.Qk(k, a, b, c, d);
            maps_k_a[std::size_t(k)][i].emplace_back(ix, qk);
          }
        }
      }
    }
  }

  std::cout << "..\n";
  std::cout << "Fill vector: " << t.reading_str() << std::endl;
  t.start();
  // std::cin.get();

  // Three-step method for filling map cut down time by >4x!

  // 1) re-size map
  m_data.resize(max_k + 1);

  // 2) Reserve enough space in each sub-map (=> 2x speed-up!)
  for (auto ik = 0ul; ik <= max_k; ++ik) {
    auto size_k = 0ul;
    for (auto ia = 0ul; ia < basis.size(); ++ia) {
      size_k += maps_k_a[ik][ia].size();
    }
    m_data[ik].reserve(size_k);
  }

// 3) Transfer data to map (can //-ize over k (since map is re-sized!))
#pragma omp parallel for
  for (auto ik = 0ul; ik <= max_k; ++ik) {
    for (const auto &k_maps : maps_k_a[ik]) {
      m_data[ik].insert(k_maps.begin(), k_maps.end());
    }
    maps_k_a[ik].clear(); // can clear vector here to "save" memory..
  }

  std::cout << "Fill map: " << t.lap_reading_str() << std::endl;
  count();
  // std::cin.get();
}

//******************************************************************************
void CoulombTable::write(const std::string &fname) const {
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
}
bool CoulombTable::read(const std::string &fname) {
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
      BigIndex key;
      Real value;
      rw_binary(f, rw, key, value);
      Q_k[key] = value;
    }
  }
  return true;
}

} // namespace Coulomb
