#pragma once
#include "Angular/Angular_369j.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace Angular {

//******************************************************************************

/*!
@brief
Lookup table for Wigner 6J symbols.

@details
Makes use of symmetry.
*/
class SixJTable {
private:
  std::vector<
      std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>
      m_data{};
  int m_max_2jk{-1};
  static auto s(int i) { return static_cast<std::size_t>(i); };

public:
  SixJTable() = default;

  //! Fill the table. max_2jk is 2* the maximum j/k that appears in the symnbols
  //! (note: 2*, as integer) - see fill()
  SixJTable(int max_2jk) { fill(max_2jk); }

  //! Returns 2* maximum k in the tales
  int max_2jk() const { return m_max_2jk; }

  //! Return 6j symbol {a/2,b/2,c/2,d/2,e/2,f/2}. Note: takes in 2*j as int
  //! @details Note: If requesting a 6J symbol beyond what is stores, will
  //! return 0 (without warning)
  double operator()(int a, int b, int c, int d, int e, int f) const {
    return get(a, b, c, d, e, f);
  }

  //! Return 6j symbol {a/2,b/2,c/2,d/2,e/2,f/2}. Note: takes in 2*j as int
  //! @details Note: If requesting a 6J symbol beyond what is stores, will
  //! return 0 (without warning)
  double get(int a, int b, int c, int d, int e, int f) const {
    if (Angular::sixj_zeroQ(a, b, c, d, e, f))
      return 0.0;
    // Faster to short-circut 0 here due to triangle? Or not?
    if (std::max({a, b, c, d, e, f}) > m_max_2jk)
      return 0.0; // error?
    const auto [i, j, k, l, m, n] = normal_order(a, b, c, d, e, f);
    return m_data[i][j][k][l][m][n];
  }

  //! Fill the table. max_2jk is 2* the maximum j/k that appears in the symnbols
  //! (note: 2*, as integer).
  /*! @details Typically, max_2jk is 2*max_2j, where max_2j is 2* max j of set
    of orbitals. You may call this function several times; if the new max_2jk is
    larger, it will extend the table. If it is smaller, does nothing. */
  void fill(int max_2jk) {
    const auto min_2jk = m_max_2jk + 1;
    // a = min{a,b,c,d,e,f}
    // b = min{b, c, e, f}

    if (max_2jk <= m_max_2jk)
      return;

    // Re-size table:
    // Note: inefficient, mostly these are empty
    // but, random access seems much faster than hash table
    const auto nsize = std::size_t(max_2jk + 1);
    m_data.resize(nsize);
    for (auto &di : m_data) {
      di.resize(nsize);
      for (auto &dij : di) {
        dij.resize(nsize);
        for (auto &dijk : dij) {
          dijk.resize(nsize);
          for (auto &dijkl : dijk) {
            dijkl.resize(nsize);
            for (auto &dijklm : dijkl) {
              dijklm.resize(nsize);
            }
          }
        }
      }
    }

    // Calculate all new *unique* 6J symbols:
    // Take advantage of symmetries, to only calc those that are needed.
    // We define a = min(a,b,c,d,e,f), b=min(b,d,e,f) => unique 6j symbol
    // in "normal" order
    for (int a = min_2jk; a <= max_2jk; ++a) {
      for (int b = a; b <= max_2jk; ++b) {
        for (int c = b; c <= max_2jk; ++c) {
          for (int d = a; d <= max_2jk; ++d) { // note: different!
            for (int e = b; e <= max_2jk; ++e) {
              for (int f = b; f <= max_2jk; ++f) {
                if (Angular::sixj_zeroQ(a, b, c, d, e, f))
                  continue;
                const auto sj = Angular::sixj_2(a, b, c, d, e, f);
                if (sj != 0.0) {
                  m_data[s(a)][s(b)][s(c)][s(d)][s(e)][s(f)] = sj;
                }
              }
            }
          }
        }
      }
    }

    // update max 2k
    m_max_2jk = max_2jk;
  }

private:
  inline static auto normal_order_level2(int a, int b, int c, int d, int e,
                                         int f) {
    // note: 'a' must be minimum!
    // assert(a == std::min({a, b, c, d, e, f})); // remove
    // {a,b,c|d,e,f} = {a,c,b|d,f,e} = {a,e,f|d,b,c} = {a,f,e|d,c,b}
    const auto min_bcef = std::min({b, c, e, f});

    if (min_bcef == b) {
      return std::array{s(a), s(b), s(c), s(d), s(e), s(f)};
    } else if (min_bcef == c) {
      return std::array{s(a), s(c), s(b), s(d), s(f), s(e)};
    } else if (min_bcef == e) {
      return std::array{s(a), s(e), s(f), s(d), s(b), s(c)};
    } else if (min_bcef == f) {
      return std::array{s(a), s(f), s(e), s(d), s(c), s(b)};
    }
    assert(false && "Fatal error 170: unreachable");
  }

  static std::array<std::size_t, 6> normal_order(int a, int b, int c, int d,
                                                 int e, int f) {
    // returns unique "normal ordering" of {a,b,c,d,e,f}->{i,j,k,l,m,n}
    // where i = min{a,b,c,d,e,f}, j = min{b,c,e,f}
    const auto min = std::min({a, b, c, d, e, f});
    //   {a,b,c|d,e,f} = {b,a,c|e,d,f} = {c,a,b|f,d,e}
    // = {d,e,c|a,b,f} = {e,d,c|b,a,f} = {f,a,e|c,d,b}
    // at next level, use also:
    // {a,b,c|d,e,f} = {a,c,b|d,f,e} = {a,e,f|d,b,c} = {a,f,e|d,c,b}
    if (min == a) {
      return normal_order_level2(a, b, c, d, e, f);
    } else if (min == b) {
      return normal_order_level2(b, a, c, e, d, f);
    } else if (min == c) {
      return normal_order_level2(c, a, b, f, d, e);
    } else if (min == d) {
      return normal_order_level2(d, e, c, a, b, f);
    } else if (min == e) {
      return normal_order_level2(e, d, c, b, a, f);
    } else if (min == f) {
      return normal_order_level2(f, a, e, c, d, b);
    }
    assert(false && "Fatal error 193: unreachable");
  }
};

} // namespace Angular
