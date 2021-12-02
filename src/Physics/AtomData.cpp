#include "AtomData.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype> //char from string
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility> // pair
#include <vector>

//******************************************************************************
namespace AtomData {

static inline bool string_is_ints(std::string_view s) {
  return !s.empty() //
         && std::find_if(s.cbegin() + 1, s.cend(),
                         [](auto c) { return !std::isdigit(c); }) == s.end() //
         && (std::isdigit(s[0]) || (s[0] == '-' && s.size() > 1));
}

//******************************************************************************
std::string int_to_roman(int a) {
  if (a > 3999)
    return std::to_string(a);
  static const std::string M[] = {"", "M", "MM", "MMM"};
  static const std::string C[] = {"",  "C",  "CC",  "CCC",  "CD",
                                  "D", "DC", "DCC", "DCCC", "CM"};
  static const std::string X[] = {"",  "X",  "XX",  "XXX",  "XL",
                                  "L", "LX", "LXX", "LXXX", "XC"};
  static const std::string I[] = {"",  "I",  "II",  "III",  "IV",
                                  "V", "VI", "VII", "VIII", "IX"};
  return M[a / 1000] + C[(a % 1000) / 100] + X[(a % 100) / 10] + I[(a % 10)];
}

//******************************************************************************
std::string NonRelSEConfig::symbol() const {
  return std::to_string(n) + AtomData::l_symbol(l) + std::to_string(num);
}

bool NonRelSEConfig::ok() const {
  if (l + 1 > n || l < 0 || n < 1 || num < 0 || num > 4 * l + 2)
    return false;
  return true;
}

//******************************************************************************
int defaultA(int Z)
// c++14: can make constexpr ?
{
  auto match_Z = [Z](const Element &atom) { return atom.Z == Z; };
  auto atom =
      std::find_if(periodic_table.begin(), periodic_table.end(), match_Z);
  if (atom == periodic_table.end())
    return 0;
  return atom->A;
}

std::string atomicSymbol(int Z) {
  auto match_Z = [Z](const Element &atom) { return atom.Z == Z; };
  auto atom =
      std::find_if(periodic_table.begin(), periodic_table.end(), match_Z);
  if (atom == periodic_table.end())
    return std::to_string(Z);
  return atom->symbol;
}

std::string atomicName(int Z) {
  auto match_Z = [Z](const Element &atom) { return atom.Z == Z; };
  auto atom =
      std::find_if(periodic_table.begin(), periodic_table.end(), match_Z);
  if (atom == periodic_table.end())
    return std::string("E") + std::to_string(Z);
  return atom->name;
}

// Given an atomic symbol (H, He, etc.), will return Z
// Note: Symbol must be exact, including capitalisation
int atomic_Z(const std::string &at) {

  if (at.empty())
    return 0;

  auto match_At = [=](const Element &atom) {
    return atom.symbol == at || '[' + atom.symbol + ']' == at;
  };
  const auto atom =
      std::find_if(periodic_table.begin(), periodic_table.end(), match_At);
  if (atom != periodic_table.end())
    return atom->Z;

  int z = 0;
  if (string_is_ints(at))
    z = std::stoi(at);
  if (z <= 0) {
    // std::cerr << "Invalid atom/Z: " << at << "\n";
    z = 0;
  }
  return z;
}

//******************************************************************************
// Short function that returns orbital term given l
std::string l_symbol(int l) {
  if (l < (int)spectroscopic_notation.length() && l >= 0)
    return spectroscopic_notation.substr(static_cast<unsigned long>(l), 1);
  else
    return "[" + std::to_string(l) + "]";
}

std::string kappa_symbol(int kappa) {
  const auto lstr = l_symbol(l_k(kappa)); //
  return lstr + "_" + std::to_string(twoj_k(kappa)) + "/2";
}

int symbol_to_l(std::string_view l_str) {
  for (auto i = 0ul; i < spectroscopic_notation.length(); i++) {
    if (spectroscopic_notation[i] == l_str[0])
      return int(i);
  }
  int l = -1;
  if (string_is_ints(l_str))
    l = std::stoi(std::string(l_str));
  else
    std::cerr << "\nFAIL AtomData::69 Invalid l: " << l_str << "?\n";

  return l;
}

//------------------------------------------------------------------------------
std::pair<int, int> parse_symbol(std::string_view symbol) {

  const auto l_ptr =
      std::find_if(symbol.begin(), symbol.end(),
                   [](const char &c) { return !std::isdigit(c); });
  const auto l_pos = std::size_t(l_ptr - symbol.begin());

  const auto n = (string_is_ints(symbol.substr(0, l_pos - 0))) ?
                     std::stoi(std::string(symbol.substr(0, l_pos - 0))) :
                     0;

  const auto l =
      (l_pos < symbol.size()) ? symbol_to_l(symbol.substr(l_pos, 1)) : -1;

  int kappa = l == 0 ? -1 : 0; // allow '6s' instead of '6s+'
  if (l >= 0 && l_pos + 1 < symbol.size()) {
    const auto pm = symbol.substr(l_pos + 1);
    const auto tj = pm == "+" ? 2 * l + 1 : pm == "-" ? 2 * l - 1 : 0;
    if (tj != 0)
      kappa = kappa_twojl(tj, l);
  }

  return {n, kappa};
}

//******************************************************************************
std::string coreConfig(const std::string &in_ng) {
  // Note: must return SAME string if no matching Nobel Gas found
  // (so that this doesn't break if I give it a full term list)
  using StringPair = std::pair<std::string, std::string>;
  auto match_ng = [&](const StringPair &ng) { return ng.first == in_ng; };
  auto ng_config =
      std::find_if(nobelGasses.begin(), nobelGasses.end(), match_ng);
  if (ng_config == nobelGasses.end())
    return in_ng;
  return ng_config->second;
}

std::string niceCoreOutput(const std::string &full_core) {
  // nb: there are 8 _actual_ nobel gasses (including H-like).
  // Only want actual nobel gasses in 'nice' output
  std::string nice_core = full_core;
  for (int i = 7; i >= 0; i--) { // loop backwards (so can break)
    auto &ng_fullterm = nobelGasses[std::size_t(i)].second;
    if (full_core.rfind(ng_fullterm, 0) == 0) {
      nice_core = nobelGasses[std::size_t(i)].first +
                  full_core.substr(ng_fullterm.length());
      break;
    }
  }

  return nice_core;
}

//******************************************************************************
double diracen(double z, double n, int k, double alpha) {
  const double a2 = alpha * alpha;
  const double c2 = 1.0 / a2;
  const double za2 = z * z * a2;
  const double g = std::sqrt(k * k - za2);
  const double w2 = z * z / std::pow(g + n - fabs((double)k), 2);
  const double d = 1.0 + a2 * w2;
  return -w2 / (2 * d) -
         (0.5 * a2 * w2 + 1. - std::sqrt(1.0 + a2 * w2)) * (c2 / d);
}

//******************************************************************************
// Takes a "core string" in form "[X],nLm,nLm,..." converts to vector of
// NonRelSEConfig, after converting [X] to a state string. Allows negative and
// non-physical m's (to allow combining); responsability of whoever uses the
// list to check for validity.
std::vector<NonRelSEConfig> core_parser(const std::string &str_core_in) {

  std::vector<NonRelSEConfig> core;
  // long unsigned int beg = 0;

  bool first = true;

  std::istringstream ss(str_core_in);
  std::string each;
  while (std::getline(ss, each, ',')) {
    if (first) {
      first = false;
      auto str_core = coreConfig(each);
      if (str_core != each) {
        state_parser(&core, str_core);
        continue;
      } else {
        auto z = atomic_Z(each);
        if (z != 0) {
          core = core_guess(z);
          continue;
        }
      }
    }
    auto new_config = term_parser(each);

    // Check if valid:
    if (new_config.n <= 0) {
      std::cout << "Problem with core: " << str_core_in << "\n";
      std::cout << "invalid core term: " << each << "\n";
      std::abort();
    }

    // If nl term already exists, add to num. Otherwise, add new term
    auto ia = std::find(core.begin(), core.end(), new_config);
    if (ia == core.end()) {
      core.push_back(new_config);
    } else {
      *ia += new_config;
    }
  }

  // std::cout << "\n----\n";
  // std::cout << str_core_in << ":\n";
  // for (auto &c : core) {
  //   std::cout << c.symbol() << ", ";
  // }
  // std::cout << "\n----\n\n";
  // std::cin.get();

  if (core.size() > 0) {
    for (auto it = core.end() - 1; it != core.begin(); it--) {
      if (it->num == 0)
        core.erase(it);
    }
  }

  return core;
}

//******************************************************************************
NonRelSEConfig term_parser(std::string_view term) {
  if (term == "" || term == "0")
    return NonRelSEConfig(0, 0, 0);

  bool term_ok = true;

  // find position of 'l'
  const auto l_ptr = std::find_if(
      term.begin(), term.end(), [](const char &c) { return !std::isdigit(c); });
  const auto l_position = std::size_t(l_ptr - term.begin());
  // Extract n, num, and l:
  int n{0}, num{-1}, l{-1};
  if (string_is_ints(term.substr(0, l_position - 0)))
    n = std::stoi(std::string(term.substr(0, l_position - 0)));
  if (term.size() > l_position + 1) {
    if (string_is_ints(term.substr(l_position + 1)))
      num = std::stoi(std::string(term.substr(l_position + 1)));
  } else {
    // string too short, mussing 'num' after l
    term_ok = false;
  }
  if (l_position == term.size())
    term_ok = false;
  if (term_ok)
    l = AtomData::symbol_to_l(term.substr(l_position, 1));

  if (num == 0)
    return NonRelSEConfig(0, 0, 0);

  // Check if valid:
  if (!term_ok || n <= 0 || l < 0) {
    return NonRelSEConfig(0, 0, 0);
  }

  // If nl term already exists, add to num. Otherwise, add new term
  return NonRelSEConfig(n, l, num);
}

//******************************************************************************
// Takes a string of states in form "nLm,nLm,..." converts to vector of
// NonRelSEConfig. Allows negative and non-physical m's (to allow combining)
std::vector<NonRelSEConfig> state_parser(const std::string &str_states) {
  std::vector<NonRelSEConfig> states;
  state_parser(&states, str_states);
  return states;
}

void state_parser(std::vector<NonRelSEConfig> *states,
                  const std::string &str_states) {

  if (str_states == "")
    return;

  std::istringstream ss(str_states);
  std::string term;
  while (std::getline(ss, term, ',')) {

    auto new_config = term_parser(term);

    // Check if valid:
    if (new_config.n <= 0) {
      std::cout << "Problem with core: " << str_states << "\n";
      std::cout << "invalid core term: " << term << "\n";
      std::abort();
    }

    // If nl term already exists, add to num. Otherwise, add new term
    auto ia = std::find(states->begin(), states->end(), new_config);
    if (ia == states->end()) {
      states->push_back(new_config);
    } else {
      *ia += new_config;
    }
  }

  // Remove any trailing states with m=zero (only if trailing; not required,
  // but nicer)
  while (!states->empty() && states->back().num == 0)
    states->pop_back();
}

//------------------------------------------------------------------------------
std::string guessCoreConfigStr(const int total_core_electrons) {

  auto core_configs = core_guess(total_core_electrons);

  static const std::vector<int> nobel_gas_list = {118, 86, 54, 36,
                                                  18,  10, 2,  0};
  static const std::vector<std::string> ng_symb_list = {
      "[Og]", "[Rn]", "[Xe]", "[Kr]", "[Ar]", "[Ne]", "[He]", "[]"};

  switch (total_core_electrons) {
    // A few over-writes for common 'different' ones
    // Lanthenides/Actinides still possibly incorrect
  case 29: // Cu
    return "[Ar],3d10,4s1";
  case 47: // Ag
    return "[Kr],4d10,5s1";
  case 79: // Au
    return "[Xe],4f14,5d10,6s1";
  }

  std::string output_config = "";
  auto index = 0u;
  for (auto n_ng : nobel_gas_list) {
    std::string ng = ng_symb_list[index++];
    if (n_ng > total_core_electrons)
      continue;
    auto ng_config = core_guess(n_ng);
    output_config += ng;
    for (auto i = ng_config.size(); i < core_configs.size(); ++i) {
      output_config += "," + core_configs[i].symbol();
    }
    break;
  }
  return output_config;
}

//------------------------------------------------------------------------------
std::vector<NonRelSEConfig> core_guess(const int total_core_electrons) {
  auto core_configs = AtomData::state_parser(AtomData::filling_order);
  auto nel = total_core_electrons;
  for (auto &c : core_configs) {
    if (c.num > nel) {
      c.num = nel;
    }
    nel -= c.num;
  }
  while (!core_configs.empty() && core_configs.back().num == 0)
    core_configs.pop_back();

  return core_configs;
}

//******************************************************************************
std::vector<DiracSEnken> listOfStates_nk(const std::string &in_list) {
  std::vector<DiracSEnken> state_list;

  std::string n_str_previous = "";
  std::string n_str = "";
  for (char c : in_list) {
    if (std::isdigit(c)) {
      n_str += c;
    } else {
      if (n_str == "")
        n_str = n_str_previous;

      const int n_max = (string_is_ints(n_str)) ? std::stoi(n_str) : -1;
      const auto l_str = std::string(1, c);
      const auto l = AtomData::symbol_to_l(l_str);

      if (l != 0) {
        for (int n = l + 1; n <= n_max; ++n)
          state_list.emplace_back(n, l);
      }
      for (int n = l + 1; n <= n_max; ++n)
        state_list.emplace_back(n, -l - 1);

      n_str_previous = n_str;
      n_str = "";
    }
  }
  return state_list;
}

//******************************************************************************
std::vector<DiracSEnken> listOfStates_singlen(const std::string &in_list) {
  std::vector<DiracSEnken> state_list;

  std::string n_str_previous = "999";
  std::string n_str = "";
  for (char c : in_list) {
    if (std::isdigit(c)) {
      n_str += c;
    } else {
      if (n_str == "")
        n_str = n_str_previous;
      const int n_max = (string_is_ints(n_str)) ? std::stoi(n_str) : -1;

      auto l_str = std::string(1, c);
      auto l = AtomData::symbol_to_l(l_str);

      if (l != 0) {
        state_list.emplace_back(n_max, l);
      }
      state_list.emplace_back(n_max, -l - 1);

      n_str_previous = n_str;
      n_str = "";
    }
  }
  return state_list;
}

//******************************************************************************
inline std::string helper_s(const Element &el) {
  auto sym = el.symbol;
  auto sym_buff = (sym.length() == 1) ? std::string("  ") : std::string(" ");
  return sym_buff + sym + " ";
}
inline std::string helper_z(const Element &el) {
  auto z_str = std::to_string(el.Z);
  auto Z_buff = (el.Z < 10) ? std::string("  ") :
                              (el.Z < 100) ? std::string(" ") : std::string("");
  return Z_buff + z_str + " ";
}

//******************************************************************************
void printTable() {

  std::string output = "";
  std::string output_s = "";
  std::string output_z = "";

  std::string output_s_l = "";
  std::string output_z_l = "";
  std::string output_s_a = "";
  std::string output_z_a = "";

  auto spaces2 = [](int n) {
    std::string buff = "";
    for (int i = 0; i < n; ++i)
      buff += "    ";
    return buff;
  };

  int row = 1;
  int col = 1;
  for (const auto &el : periodic_table) {
    if (el.Z > 118)
      break;

    if (row == 1 && col == 2) {
      col += 16;
      output_s += spaces2(16);
      output_z += spaces2(16);
    } else if (row < 4 && col == 3) {
      col += 10;
      output_s += spaces2(10);
      output_z += spaces2(10);
    } else if (row > 5) {
      if (col == 2) {
        output_s += " *  ";
        output_z += "    ";
        col++;
      }
      if (el.Z >= 57 && el.Z < 72) {
        output_s_l += helper_s(el);
        output_z_l += helper_z(el);
        continue;
      } else if (el.Z >= 89 && el.Z < 104) {
        output_s_a += helper_s(el);
        output_z_a += helper_z(el);
        continue;
      }
    }

    output_s += helper_s(el);
    output_z += helper_z(el);
    ++col;

    if (col > 17 || el.Z == 118) {
      ++row;
      col = 0;
      // output += std::string(output_s + "\n" + output_z + "\n");
      output += output_s;
      output += "\n";
      output += output_z;
      output += "\n";
      output_s.clear();
      output_z.clear();
    }
  }
  std::cout << "\n" << output << "\n";
  std::cout << "      * " << output_s_l << "\n";
  std::cout << "        " << output_z_l << "\n";
  std::cout << "      * " << output_s_a << "\n";
  std::cout << "        " << output_z_a << "\n\n";
}

} // namespace AtomData
