#pragma once
#include "Angular/Wigner369j.hpp"
#include "Wavefunction/DiracSpinor.hpp" // for 'magic' 6J symbols
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include <vector>

// XXX Note: This is significantly faster if implemented in header file, not
// sepperate cpp file. Seems due to inlineing of 'get' function

namespace Angular {

//******************************************************************************

/*!
@brief
Lookup table for Wigner 6J symbols.

@details
Note: functions all called with 2*j and 2*k (ensure j integer)
Makes use of symmetry.
*/
class SixJTable {
private:
  std::unordered_map<uint64_t, double> m_data{};
  int m_max_2j_k{-1};
  static auto s(int i) { return static_cast<uint8_t>(i); };

public:
  SixJTable() = default;

  //! Fill the table. max_2j_k is the maximum k that appears in the symbols,
  //! equivilant to 2*max(j) (always integer) - see fill()
  SixJTable(int max_2j_k) { fill(max_2j_k); }

  //! Returns maximum k in the tables [max(k)=2*max(j)=max(2j)]
  int max_2jk() const { return m_max_2j_k; }

  //----------------------------------------------------------------------------
  //! Return 6j symbol {a/2,b/2,c/2,d/2,e/2,f/2}. Note: takes in 2*j/2*k as int
  //! @details Note: If requesting a 6J symbol beyond what is stored, will
  //! return 0 (without warning)
  inline double get_2(int a, int b, int c, int d, int e, int f) const {
    if (Angular::sixj_zeroQ(a, b, c, d, e, f))
      return 0.0;
    const auto it = m_data.find(normal_order(a, b, c, d, e, f));
    return (it == m_data.cend()) ? 0.0 : it->second;
  }

  //! Return "magic" 6j. Pass in either integer (for k), or DiracSpinor, F.
  //! e.g.: (F,F,k,F,F,k) -> {j,j,k,j,j,k}. Do NOT *2!
  //! @details Note: If requesting a 6J symbol beyond what is stored, will
  //! return 0 (without warning)
  template <class A, class B, class C, class D, class E, class F>
  double get(const A &a, const B &b, const C &c, const D &d, const E &e,
             const F &f) const {
    return get_2(twojk(a), twojk(b), twojk(c), twojk(d), twojk(e), twojk(f));
  }

  //----------------------------------------------------------------------------
  //! Checks if given 6j symbol is in table (note: may not be in table because
  //! it's zero)
  bool contains(int a, int b, int c, int d, int e, int f) const {
    const auto it = m_data.find(normal_order(a, b, c, d, e, f));
    return (it != m_data.cend());
  }

  //----------------------------------------------------------------------------
  //! Fill the table. max_2jk is 2* the maximum j/k that appears in the symnbols
  //! (note: 2*, as integer).
  /*! @details Typically, max_2jk is 2*max_2j, where max_2j is 2* max j of set
    of orbitals. You may call this function several times; if the new max_2jk is
    larger, it will extend the table. If it is smaller, does nothing. */
  void fill(int max_2j_k) {
    // a = min{a,b,c,d,e,f}
    // b = min{b, c, e, f}

    if (max_2j_k <= m_max_2j_k)
      return;

    // Calculate all new *unique* 6J symbols:
    // Take advantage of symmetries, to only calc those that are needed.
    // We define a = min(a,b,c,d,e,f), b=min(b,d,e,f) => unique 6j symbol
    // in "normal" order
    const auto max_2k = 2 * max_2j_k;
    for (int a = 0; a <= max_2k; ++a) {
      auto a0 = a; // std::max(min_2jk, a);
      for (int b = a0; b <= max_2k; ++b) {
        for (int c = b; c <= max_2k; ++c) {
          for (int d = a0; d <= max_2k; ++d) { // note: different!
            for (int e = b; e <= max_2k; ++e) {
              for (int f = b; f <= max_2k; ++f) {
                if (Angular::sixj_zeroQ(a, b, c, d, e, f))
                  continue;
                if (contains(a, b, c, d, e, f))
                  continue;
                const auto sj = Angular::sixj_2(a, b, c, d, e, f);
                if (std::abs(sj) > 1.0e-16) {
                  m_data[normal_order(a, b, c, d, e, f)] = sj;
                }
              }
            }
          }
        }
      }
    }

    // update max 2k
    m_max_2j_k = max_2j_k;
  }

private:
  //----------------------------------------------------------------------------
  inline static auto make_key(uint8_t a, uint8_t b, uint8_t c, uint8_t d,
                              uint8_t e, uint8_t f) {
    uint64_t key = 0;
    // std::memcpy(&key, &abcdef[0], 6 * sizeof(uint8_t));
    auto pk = (uint8_t *)(&key);
    std::memcpy(pk, &a, sizeof(uint8_t));
    std::memcpy(pk + 1, &b, sizeof(uint8_t));
    std::memcpy(pk + 2, &c, sizeof(uint8_t));
    std::memcpy(pk + 3, &d, sizeof(uint8_t));
    std::memcpy(pk + 4, &e, sizeof(uint8_t));
    std::memcpy(pk + 5, &f, sizeof(uint8_t));
    return key;
  }

  //----------------------------------------------------------------------------
  inline static auto normal_order_level2(int a, int b, int c, int d, int e,
                                         int f) {
    // note: 'a' must be minimum!
    // assert(a == std::min({a, b, c, d, e, f})); // remove
    // {a,b,c|d,e,f} = {a,c,b|d,f,e} = {a,e,f|d,b,c} = {a,f,e|d,c,b}
    const auto min_bcef = std::min({b, c, e, f});

    if (min_bcef == b) {
      return make_key(s(a), s(b), s(c), s(d), s(e), s(f));
    } else if (min_bcef == c) {
      return make_key(s(a), s(c), s(b), s(d), s(f), s(e));
    } else if (min_bcef == e) {
      return make_key(s(a), s(e), s(f), s(d), s(b), s(c));
    } else if (min_bcef == f) {
      return make_key(s(a), s(f), s(e), s(d), s(c), s(b));
    }
    assert(false && "Fatal error 170: unreachable");
  }

  //----------------------------------------------------------------------------
  static uint64_t normal_order(int a, int b, int c, int d, int e, int f) {
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

  //******************************************************************************
  // If given an integer (k), returns 2*k
  // If given a DiracSpinor, F, returns 2*j [F.twoj()]
  template <class A> static int twojk(const A &a) {
    if constexpr (std::is_same_v<A, DiracSpinor>) {
      return a.twoj();
    } else {
      static_assert(std::is_same_v<A, int>);
      return 2 * a;
    }
  }
};

} // namespace Angular
