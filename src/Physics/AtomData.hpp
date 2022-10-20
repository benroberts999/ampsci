#pragma once
#include "Angular/Wigner369j.hpp"
#include "Physics/AtomData_PeriodicTable.hpp"
#include "qip/Template.hpp"
#include <string>
#include <utility>
#include <vector>

//! Useful atomic data/functions
namespace AtomData {

//==============================================================================
//! Stores non-relativistic single-eletron config: {n, l, number}
struct NonRelConfig : qip::Comparison<NonRelConfig>,
                      qip::Arithmetic<NonRelConfig> {
  int n;
  int l;
  int num;

  constexpr NonRelConfig(int in_n = 0, int in_l = -1, int in_num = 0)
      : n(in_n), l(in_l), num(in_num) {}

  //! Returns symbol (e.g., 1s2 or 5p3)
  std::string symbol() const;

  //! Checks if consistent (l>n etc.)
  bool ok() const;

  //! Filling fraction (accounting for spin) = num/[2*(2l+1)]
  double frac() const;

  //! Provides comparitor overloads. Compares n first, then l.
  friend bool operator==(const NonRelConfig &lhs, const NonRelConfig &rhs);
  friend bool operator<(const NonRelConfig &lhs, const NonRelConfig &rhs);

  //! Provides addition and subtraction: adds 'num' iff n and l same
  NonRelConfig &operator+=(const NonRelConfig &rhs);
  NonRelConfig &operator-=(const NonRelConfig &rhs);
};

//==============================================================================
//! Stores relativistic single-eletron state {n, kappa, energy}
struct DiracConfig { // name OK? too short?
  int n;
  int k;
  double en;
  DiracConfig(int in_n = 0, int in_k = 0, double in_en = 0)
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

//!  Returns shortSymbol, given n and kappa: (6,-1)->"6s+"
std::string shortSymbol(int n, int kappa);

//! Parses electron 'symbol' or 'shortSymbol' to {n,kappa}, e.g., "6s+" -> {6,-1}; "6p-" -> {6,1}; "6p_1/2" -> {6,1}
std::pair<int, int> parse_symbol(std::string_view symbol);

//! Exact H-like energy
double diracen(double z, double n, int k, double alpha = 0.00729735256635);

//! Prints a periodic table to screen
void printTable();

//! Given a nobel-gas conifg (e.g., '[Xe]') returns full electron config
std::string coreConfig(const std::string &in_ng);

//! Given a full electron config., returns nicer format by recognising nobel gas
std::string niceCoreOutput(const std::string &full_core);

//! Takes a "core string" in form "[X],nLm,nLm,..." converts to vector of
//! NonRelConfig, after converting [X] to a state string. Allows negative and
//! non-physical m's (to allow combining); responsability of whoever uses the
//! list to check for validity.
std::vector<NonRelConfig> core_parser(const std::string &str_core_in);

//! Given a list of NonRelConfigs, returns full string
std::string configs_to_string(const std::vector<NonRelConfig> &configs);

//! Given a term symbol 'nLm', returns corresponding NonRelConfig
NonRelConfig term_parser(std::string_view term);

//! Takes a string of states in form "nLm,nLm,..." converts to vector of
//! NonRelConfig. Allows negative and non-physical m's (to allow combining)
std::vector<NonRelConfig> state_parser(const std::string &str_states);
//! Overload; adds to existing states vector (may be empty)
void state_parser(std::vector<NonRelConfig> *states,
                  const std::string &str_states);

//! Given a number of electrons, guesses the configuration, returns as string
std::string guessCoreConfigStr(const int total_core_electrons);

//! Given a number of electrons, guesses the configuration, returns list of NonRelConfigs
std::vector<NonRelConfig> core_guess(const int total_core_electrons);

//! Generates a list of DiracConfig from string: full list
std::vector<DiracConfig> listOfStates_nk(const std::string &in_list);

//! Generates a list of DiracConfig from string: just max n for each kappa
std::vector<DiracConfig> listOfStates_singlen(const std::string &in_list);

} // namespace AtomData
