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
void CoulombTable::summary() const {
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

int CoulombTable::count() const {
  auto total = 0ul;
  for (auto &qk : m_data) {
    total += qk.size();
  }
  return static_cast<int>(total);
}

//******************************************************************************
void CoulombTable::add(int k, const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d, Real value) {
  // [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // const auto sk = std::size_t(k);
  // if (sk >= m_data.size()) {
  //   m_data.resize(sk + 1);
  // }
  // m_data.at(sk).insert({NormalOrder(a, b, c, d), value});
  add(k, NormalOrder(a, b, c, d), value);
}
//----------------
void CoulombTable::add(int k, BigIndex index, Real value) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__, "2");
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
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // const auto sk = std::size_t(k);
  // if (sk >= m_data.size()) {
  //   m_data.resize(sk + 1);
  // }
  // m_data.at(sk).insert_or_assign(NormalOrder(a, b, c, d), value);
  update(k, NormalOrder(a, b, c, d), value);
}

void CoulombTable::update(int k, BigIndex index, Real value) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  const auto sk = std::size_t(k);
  if (sk >= m_data.size()) {
    m_data.resize(sk + 1);
  }
  m_data.at(sk).insert_or_assign(index, value);
}
//******************************************************************************
bool CoulombTable::contains(int k, const DiracSpinor &a, const DiracSpinor &b,
                            const DiracSpinor &c, const DiracSpinor &d) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return false;
  return m_data.at(sk).find(NormalOrder(a, b, c, d)) != m_data.at(sk).cend();
}

bool CoulombTable::contains(int k, BigIndex index) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return false;
  return m_data.at(sk).find(index) != m_data.at(sk).cend();
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
//
double CoulombTable::Q(int k, BigIndex index) const {
  const auto sk = std::size_t(k);
  if (sk >= m_data.size())
    return 0.0;
  // check valid k? Probably faster to lookup in table
  const auto map_it = m_data.at(sk).find(index);
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
                       const DiracSpinor &c, const DiracSpinor &d,
                       const Angular::SixJTable *const sj) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  double Pk_abcd{0.0};

  // 6j(s) Triads: {a,c,k}, {k,b,d}, {c,b,l}, {d,a,l}
  if (Angular::triangle(a.twoj(), c.twoj(), 2 * k) == 0)
    return 0.0;
  if (Angular::triangle(2 * k, b.twoj(), d.twoj()) == 0)
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

//******************************************************************************
double CoulombTable::W(int k, const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d,
                       const Angular::SixJTable *const sj) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  return Q(k, a, b, c, d) + P(k, a, b, c, d, sj);
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
CoulombTable::BigIndex CoulombTable::NormalOrder(const DiracSpinor &a,
                                                 const DiracSpinor &b,
                                                 const DiracSpinor &c,
                                                 const DiracSpinor &d) const {
  return NormalOrder_impl(a.nk_index(), b.nk_index(), c.nk_index(),
                          d.nk_index());
}

//******************************************************************************
CoulombTable::BigIndex CoulombTable::CurrentOrder(const DiracSpinor &a,
                                                  const DiracSpinor &b,
                                                  const DiracSpinor &c,
                                                  const DiracSpinor &d) const {
  return FormIndex(a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
}

//******************************************************************************
CoulombTable::BigIndex QkTable::NormalOrder_impl(Index a, Index b, Index c,
                                                 Index d) const {
  // put smallest first
  const auto min = std::min({a, b, c, d});
  if (min == a) {
    // options are abcd, and adcb
    return (b < d) ? FormIndex(a, b, c, d) : FormIndex(a, d, c, b);
  } else if (min == b) {
    // options are badc, and bcda
    return (a < c) ? FormIndex(b, a, d, c) : FormIndex(b, c, d, a);
  } else if (min == c) {
    // options are cbad, and cdab
    return (b < d) ? FormIndex(c, b, a, d) : FormIndex(c, d, a, b);
  } else if (min == d) {
    // options are dabc, and dcba
    return (a < c) ? FormIndex(d, a, b, c) : FormIndex(d, c, b, a);
  }
  assert(false);
}
//------------------------------------------------------------------------------
CoulombTable::BigIndex WkTable::NormalOrder_impl(Index a, Index b, Index c,
                                                 Index d) const {
  // put smallest first
  const auto min = std::min({a, b, c, d});
  if (min == a) {
    return FormIndex(a, b, c, d);
  } else if (min == b) {
    return FormIndex(b, a, d, c);
  } else if (min == c) {
    return FormIndex(c, d, a, b);
  } else if (min == d) {
    return FormIndex(d, c, b, a);
  }
  assert(false);
}

//------------------------------------------------------------------------------
CoulombTable::BigIndex LkTable::NormalOrder_impl(Index a, Index b, Index c,
                                                 Index d) const {
  // [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  if (a <= b) { // a = std::min(a, b);
    return FormIndex(a, b, c, d);
  } else {
    return FormIndex(b, a, d, c);
  }
  assert(false);
}

//------------------------------------------------------------------------------
CoulombTable::BigIndex NkTable::NormalOrder_impl(Index a, Index b, Index c,
                                                 Index d) const {
  return FormIndex(a, b, c, d);
}

//******************************************************************************
CoulombTable::BigIndex CoulombTable::FormIndex(Index a, Index b, Index c,
                                               Index d) const { //
  static_assert(sizeof(BigIndex) == 4 * sizeof(Index));
  BigIndex big_index;
  auto p_index = (Index *)(&big_index);
  std::memcpy(p_index, &a, sizeof(Index));
  std::memcpy(p_index + 1, &b, sizeof(Index));
  std::memcpy(p_index + 2, &c, sizeof(Index));
  std::memcpy(p_index + 3, &d, sizeof(Index));
  return big_index;
}

//******************************************************************************
std::array<CoulombTable::Index, 4>
CoulombTable::UnFormIndex(const BigIndex &index) const {
  static_assert(sizeof(BigIndex) == sizeof(std::array<Index, 4>));
  std::array<Index, 4> set;
  std::memcpy(&set, &index, sizeof(BigIndex));
  return set;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
void QkTable::fill_old(const std::vector<DiracSpinor> &basis,
                       const YkTable &yk) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  IO::ChronoTimer t("fill");

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
  for (auto ia = 0ul; ia < basis.size(); ++ia) {
    const auto &a = basis[ia];
    for (const auto &b : basis) {
      for (const auto &c : basis) {
        for (const auto &d : basis) {

          // enfore symmetry here, to avoid calculating anything twice
          const auto ix = NormalOrder(a, b, c, d);
          if (ix != CurrentOrder(a, b, c, d))
            continue;

          const auto [kmin, kmax] = k_minmax_Q(a, b, c, d);
          for (int k = kmin; k <= kmax; k += 2) {
            const auto qk = yk.Q(k, a, b, c, d);
            maps_k_a[std::size_t(k)][ia].emplace_back(ix, qk);
          }
        }
      }
    }
  }

  std::cout << "..\n";
  std::cout << "Fill vector: " << t.reading_str() << std::endl;
  t.start();

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
    for (auto &k_maps : maps_k_a[ik]) {
      m_data[ik].insert(k_maps.begin(), k_maps.end());
    }
    // maps_k_a[ik].clear(); // can clear vector here to "save" memory..
  }

  std::cout << "Fill map: " << t.lap_reading_str() << std::endl;
  summary();
}

//*****************************
void QkTable::fill(const std::vector<DiracSpinor> &basis, const YkTable &yk) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
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

  const auto max_k = std::size_t(DiracSpinor::max_tj(basis)) - 1;
  m_data.resize(max_k + 1);

  // 1) Count non-zero Q integrals (each k). Use this to 'reserve' map space
  t.start();
  std::vector<std::size_t> count_non_zero_k(max_k + 1);
#pragma omp parallel for
  for (auto k = 0; k <= int(max_k); ++k) {
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        if (b < a) // symmetry
          continue;
        for (const auto &c : basis) {
          if (c < a) // symmetry
            continue;
          if (!Angular::Ck_kk_SR(k, a.k, c.k))
            continue;
          for (const auto &d : basis) {
            if (d < b) // symmetry
              continue;
            if (!Angular::Ck_kk_SR(k, b.k, d.k))
              continue;
            ++count_non_zero_k[std::size_t(k)];
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
  for (auto k = 0; k <= int(max_k); ++k) {
    for (const auto &a : basis) {
      for (const auto &b : basis) {
        if (b < a) // symmetry
          continue;
        for (const auto &c : basis) {
          if (c < a) // symmetry
            continue;
          if (!Angular::Ck_kk_SR(k, a.k, c.k))
            continue;
          for (const auto &d : basis) {
            if (d < b) // symmetry
              continue;
            if (!Angular::Ck_kk_SR(k, b.k, d.k))
              continue;
            add(k, a, b, c, d, 0.0);
          }
        }
      }
    }
  }
  std::cout << "Fill w/ zeros: " << t.lap_reading_str() << std::endl;

  // 4) Fill the pre-constructed map with values, in parallel. Since we are not
  // adding any new elements to map, and since we are guarenteed to only access
  // each map element once, we can do this part in parallel.
  // nb: This //isation is not very efficient, though in theory it can be 100%
  t.start();
#pragma omp parallel for
  for (auto ia = 0ul; ia < basis.size(); ++ia) {
    const auto &a = basis[ia];
    for (const auto &b : basis) {
      if (b < a) // symmetry
        continue;
      for (const auto &c : basis) {
        if (c < a) // symmetry
          continue;
        for (const auto &d : basis) {
          if (d < b) // symmetry
            continue;
          const auto [kmin, kmax] = k_minmax_Q(a, b, c, d);
          for (int k = kmin; k <= kmax; k += 2) {
            update(k, a, b, c, d, yk.Q(k, a, b, c, d));
          }
        }
      }
    }
  }

  std::cout << "Fill w/ values: " << t.lap_reading_str() << std::endl;

  summary();
}

//******************************************************************************
void CoulombTable::write(const std::string &fname) const {
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
bool CoulombTable::read(const std::string &fname) {
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
      BigIndex key;
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
