#pragma once
#include "Angular/Wigner369j.hpp"
#include "Physics/AtomData_PeriodicTable.hpp"
#include <array>
#include <cmath>
#include <string>
#include <vector>

//! Useful atomic data/functions. Most self-explanatory
namespace AtomData {

//==============================================================================
//! Stores none relativistic single-eletron config {n, l, number}
struct NonRelSEConfig {
  int n;
  int l;
  int num;
  NonRelSEConfig(int in_n = 0, int in_l = -1, int in_num = 0)
      : n(in_n), l(in_l), num(in_num) {}

  std::string symbol() const;
  bool ok() const;

  double frac() const {
    int filling = 2 * (2 * l + 1);
    return (num < filling) ? double(num) / double(filling) : 1.0;
  };

  // comparitor overloads:
  bool operator==(const NonRelSEConfig &other) const {
    return n == other.n && l == other.l;
  }
  bool operator!=(const NonRelSEConfig &other) const {
    return !(*this == other);
  }
  NonRelSEConfig &operator+=(const NonRelSEConfig &other) {
    this->num += other.num;
    return *this;
  }
};

//==============================================================================
//! Stores relativistic single-eletron state {n, kappa, energy}
struct DiracSEnken { // name OK? too short?
  int n;
  int k;
  double en;
  DiracSEnken(int in_n = 0, int in_k = 0, double in_en = 0)
      : n(in_n), k(in_k), en(in_en){};
};

//==============================================================================
//==============================================================================

//! Looks up default A (most common) for given Z
int defaultA(int Z);

//! e.g., 55 -> "Cs"
std::string atomicSymbol(int Z);
//! e.g., 55 -> "Cesium"
std::string atomicName(int Z);

//! Converts atomic symbol to integer Z (e.g., 'Cs' to 55 )
int atomic_Z(const std::string &at);
//! Overload, to can call with int anyway
inline int atomic_Z(int z) { return z; }

//! l (int) to symbol (e.g., 0->'s', 1->'p')
std::string l_symbol(int l);
//! e.g., 'p' -> 1
int symbol_to_l(std::string_view l_str);
//! kappa (int) to symbol, e.g., -1 -> s_1/2
std::string kappa_symbol(int kappa);

//! Parses "short symbol" to {n,kappa}, e.g., "6s+" -> {6,-1}; "6p-" -> {6,1}
std::pair<int, int> parse_symbol(std::string_view symbol);

//! Given a nobel-gas conifg (e.g., '[Xe]') returns full electron config
std::string coreConfig(const std::string &in_ng);

//! Given a full electron config., returns nicer format by recognising nobel gas
std::string niceCoreOutput(const std::string &full_core);

//! Exact H-like energy
double diracen(double z, double n, int k, double alpha = 0.00729735256635);

//! Takes a "core string" in form "[X],nLm,nLm,..." converts to vector of
//! NonRelSEConfig, after converting [X] to a state string. Allows negative and
//! non-physical m's (to allow combining); responsability of whoever uses the
//! list to check for validity.
std::vector<NonRelSEConfig> core_parser(const std::string &str_core_in);

NonRelSEConfig term_parser(std::string_view term);

//! Takes a string of states in form "nLm,nLm,..." converts to vector of
//! NonRelSEConfig. Allows negative and non-physical m's (to allow combining)
std::vector<NonRelSEConfig> state_parser(const std::string &str_states);
//! Overload; adds to existing states vector (may be empty)
void state_parser(std::vector<NonRelSEConfig> *states,
                  const std::string &str_states);

std::string guessCoreConfigStr(const int total_core_electrons);
std::vector<NonRelSEConfig> core_guess(const int total_core_electrons);

//! Generates a list of DiracSEnken from string: full list
std::vector<DiracSEnken> listOfStates_nk(const std::string &in_list);
//! Generates a list of DiracSEnken from string: just max n for each kappa
std::vector<DiracSEnken> listOfStates_singlen(const std::string &in_list);

//! Prints a periodic table to screen
void printTable();

//! converts into to lc romain numerals
std::string int_to_roman(int a);

//==============================================================================
//! Returns number of possible states _below_ given n
constexpr int states_below_n(int n) { return n * n - 2 * n + 1; }

//! return nk_index given {n, kappa}: nk_index(n,k) := n^2 - 2n + 1 +
//! kappa_index
/*! @details   nk_index:
 For easy array access, define 1-to-1 index for each {n, kappa}:
 nk_index(n,k) := n^2 - 2n + 1 + kappa_index.
 nb: n^2 - 2n + 1 = states_below_n - number of possible states with n'<n.
 Note: ONLY valid for n >= 1 (i.e., cannot be used for general basis states)
*/
constexpr int nk_to_index(int n, int k) {
  return states_below_n(n) + Angular::indexFromKappa(k);
}

//! return {n, kappa} given nk_index:
inline std::pair<int, int> index_to_nk(int index) {
  // Better way? isqrt?
  const auto n = 1 + int(std::sqrt(index + 0.01));
  // int n = 1 + int_sqrt(index);
  const auto kappa_index = index - states_below_n(n);
  return {n, Angular::kappaFromIndex(kappa_index)};
}

} // namespace AtomData
