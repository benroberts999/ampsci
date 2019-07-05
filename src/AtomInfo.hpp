#pragma once
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace AtomInfo {

//******************************************************************************
struct Element {
  int Z;
  std::string symbol;
  int A;
  Element(int inZ, std::string insymbol, int inA)
      : Z(inZ), symbol(insymbol), A(inA) {}
};

// Default values for A for each atom.
// Goes up to E120 (Z=120)
static const std::vector<Element> periodic_table = {
    {1, "H", 1},      {2, "He", 4},     {3, "Li", 7},       {4, "Be", 9},
    {5, "B", 11},     {6, "C", 12},     {7, "N", 14},       {8, "O", 16},
    {9, "F", 19},     {10, "Ne", 20},   {11, "Na", 23},     {12, "Mg", 24},
    {13, "Al", 27},   {14, "Si", 28},   {15, "P", 31},      {16, "S", 32},
    {17, "Cl", 35},   {18, "Ar", 40},   {19, "K", 39},      {20, "Ca", 40},
    {21, "Sc", 45},   {22, "Ti", 48},   {23, "V", 51},      {24, "Cr", 52},
    {25, "Mn", 55},   {26, "Fe", 56},   {27, "Co", 59},     {28, "Ni", 59},
    {29, "Cu", 64},   {30, "Zn", 65},   {31, "Ga", 70},     {32, "Ge", 73},
    {33, "As", 75},   {34, "Se", 79},   {35, "Br", 80},     {36, "Kr", 84},
    {37, "Rb", 85},   {38, "Sr", 88},   {39, "Y", 89},      {40, "Zr", 91},
    {41, "Nb", 93},   {42, "Mo", 96},   {43, "Tc", 97},     {44, "Ru", 101},
    {45, "Rh", 103},  {46, "Pd", 106},  {47, "Ag", 108},    {48, "Cd", 112},
    {49, "In", 115},  {50, "Sn", 119},  {51, "Sb", 122},    {52, "Te", 128},
    {53, "I", 127},   {54, "Xe", 131},  {55, "Cs", 133},    {56, "Ba", 137},
    {57, "La", 139},  {58, "Ce", 140},  {59, "Pr", 141},    {60, "Nd", 144},
    {61, "Pm", 145},  {62, "Sm", 150},  {63, "Eu", 152},    {64, "Gd", 157},
    {65, "Tb", 159},  {66, "Dy", 162},  {67, "Ho", 165},    {68, "Er", 167},
    {69, "Tm", 169},  {70, "Yb", 173},  {71, "Lu", 175},    {72, "Hf", 178},
    {73, "Ta", 181},  {74, "W", 184},   {75, "Re", 186},    {76, "Os", 190},
    {77, "Ir", 192},  {78, "Pt", 195},  {79, "Au", 197},    {80, "Hg", 201},
    {81, "Tl", 204},  {82, "Pb", 207},  {83, "Bi", 209},    {84, "Po", 209},
    {85, "At", 210},  {86, "Rn", 222},  {87, "Fr", 223},    {88, "Ra", 226},
    {89, "Ac", 227},  {90, "Th", 232},  {91, "Pa", 231},    {92, "U", 238},
    {93, "Np", 237},  {94, "Pu", 244},  {95, "Am", 243},    {96, "Cm", 247},
    {97, "Bk", 247},  {98, "Cf", 251},  {99, "Es", 252},    {100, "Fm", 257},
    {101, "Md", 258}, {102, "No", 259}, {103, "Lr", 262},   {104, "Rf", 267},
    {105, "Db", 270}, {106, "Sg", 269}, {107, "Bh", 270},   {108, "Hs", 270},
    {109, "Mt", 278}, {110, "Ds", 281}, {111, "Rg", 281},   {112, "Cn", 285},
    {113, "Nh", 286}, {114, "Fl", 289}, {115, "Mc", 289},   {116, "Lv", 293},
    {117, "Ts", 293}, {118, "Og", 294}, {119, "E119", 315}, {120, "E120", 320}};

inline int defaultA(int Z)
// c++14: can make constexpr ?
{
  for (auto &atom : periodic_table) {
    if (atom.Z == Z)
      return atom.A;
  }
  return 0;
}

inline std::string atomicSymbol(int Z) {
  for (auto &atom : periodic_table) {
    if (atom.Z == Z)
      return atom.symbol;
  }
  return std::to_string(Z);
}

// Given an atomic symbol (H, He, etc.), will return Z
// Note: Symbol must be exact, including capitalisation
inline int get_z(const std::string &at) {
  for (auto &atom : periodic_table) {
    if (atom.symbol == at)
      return atom.Z;
  }
  int z = 0;
  try {
    z = std::stoi(at);
  } catch (...) {
  }
  if (z <= 0) {
    std::cerr << "Invalid atom/Z: " << at << "\n";
    std::abort();
  }
  return z;
}

//******************************************************************************
static const std::string spectroscopic_notation = "spdfghiklmnoqrtuvwxyzabc";
static const std::string Spectroscopic_Notation = "SPDFGHIKLMNOQRTUVWXYZABC";

// Short function that returns orbital term given l
inline std::string l_symbol(int l) {
  if (l < (int)spectroscopic_notation.length() && l >= 0)
    return spectroscopic_notation.substr(l, 1);
  else
    return "[" + std::to_string(l) + "]";
}

inline int symbol_to_l(const std::string &l_str) {
  for (int i = 0; i < (int)spectroscopic_notation.length(); i++) {
    if (spectroscopic_notation.substr(i, 1) == l_str)
      return i;
  }
  int l = -1;
  try {
    // Can work if given an int as a string:
    l = std::stoi(l_str);
  } catch (...) { // don't abort here (might get nice error message later)
    std::cerr << "\nFAIL AtomInfo::69 Invalid l: " << l_str << "?\n";
  }
  return l;
}

//******************************************************************************
constexpr int l_k(int ka) { return (ka > 0) ? ka : -ka - 1; }
constexpr int twoj_k(int ka) { return (ka > 0) ? 2 * ka - 1 : -2 * ka - 1; }
constexpr double j_k(int ka) {
  return (ka > 0) ? double(ka) - 0.5 : double(-ka) - 0.5;
}
constexpr int parity_k(int ka) {
  return (ka % 2 == 0) ? ((ka > 0) ? 1 : -1) : ((ka < 0) ? 1 : -1);
}
constexpr int l_tilde_k(int ka) {
  // "Complimentary l (l for lower component)"
  // l-tilde = (2j-l) = l +/- 1, for j = l +/- 1/2
  return (ka > 0) ? ka - 1 : -ka;
}
constexpr int kappa_twojl(int twoj, int l) {
  return ((2 * l - twoj) * (twoj + 1)) / 2;
}
//******************************************************************************
//    Kappa Index:
// For easy array access, define 1-to-1 index for each kappa:
// kappa: -1  1 -2  2 -3  3 -4  4 ...
// index:  0  1  2  3  4  5  6  7 ...
// kappa(i) = (-1,i+1)*(int(i/2)+1)
constexpr int indexFromKappa(int ka) {
  return (ka < 0) ? -2 * ka - 2 : 2 * ka - 1;
}
constexpr int kappaFromIndex(int i) {
  return (i % 2 == 0) ? -(i + 2) / 2 : (i + 1) / 2;
}
constexpr int twojFromIndex(int i) { return (i % 2 == 0) ? i + 1 : i; }
constexpr int lFromIndex(int i) { return (i % 2 == 0) ? i / 2 : (i + 1) / 2; }
//******************************************************************************

// Note: this requires that all Nobel Gasses are listed FIRST, in order
// (Assumed by "niceCoreOutput" function that this matches nobelGasses
static const std::array<std::pair<std::string, std::string>, 12> nobelGasses = {
    std::make_pair("[He]", "1s2"),
    /**/ //
    std::make_pair("[Ne]", "1s2,2s2,2p6"),
    std::make_pair("[Ar]", "1s2,2s2,2p6,3s2,3p6"),
    std::make_pair("[Kr]", "1s2,2s2,2p6,3s2,3p6,3d10,4s2,4p6"),
    std::make_pair("[Xe]", "1s2,2s2,2p6,3s2,3p6,3d10,4s2,4p6,4d10,5s2,5p6"),
    std::make_pair(
        "[Rn]",
        "1s2,2s2,2p6,3s2,3p6,3d10,4s2,4p6,4d10,5s2,5p6,4f14,5d10,6s2,6p6"),
    std::make_pair("[Og]", "1s2,2s2,2p6,3s2,3p6,3d10,4s2,4p6,4d10,5s2,5p6,"
                           "4f14,5d10,6s2,6p6,5f14,6d10,7s2,7p6"),
    // A few extra:
    std::make_pair("[Zn]", "1s2,2s2,2p6,3s2,3p6,3d10,4s2"),
    std::make_pair("[Cd]", "1s2,2s2,2p6,3s2,3p6,3d10,4s2,4p6,4d10,5s2"),
    std::make_pair(
        "[Hg]", "1s2,2s2,2p6,3s2,3p6,3d10,4s2,4p6,4d10,5s2,5p6,4f14,5d10,6s2"),
    std::make_pair("[Cn]", "1s2,2s2,2p6,3s2,3p6,3d10,4s2,4p6,4d10,5s2,5p6,"
                           "4f14,5d10,6s2,6p6,5f14,6d10,7s2"),
    std::make_pair("[]", "1s0")};

inline std::string coreConfig(const std::string &in_ng) {
  // Note: must return SAME string if no matching Nobel Gas found
  // (so that this doesn't break if I give it a full term list)
  for (auto &ng : nobelGasses) {
    if (in_ng == ng.first)
      return ng.second;
  }
  return in_ng;
}

inline std::string niceCoreOutput(const std::string &full_core) {
  // nb: there are 7 _actual_ nobel gasses.
  // Only want actual nobel gasses in 'nice' output
  std::string nice_core = full_core;
  for (int i = 6; i >= 0; i--) { // loop backwards (so can break)
    auto &ng_fullterm = nobelGasses[i].second;
    if (full_core.rfind(ng_fullterm, 0) == 0) {
      nice_core = nobelGasses[i].first + full_core.substr(ng_fullterm.length());
      break;
    }
  }
  return nice_core;
}

//******************************************************************************
inline double diracen(double z, double n, int k,
                      double alpha = 0.00729735256635) {
  double a2 = alpha * alpha;
  double c2 = 1. / a2;
  double za2 = z * z * a2;
  double g = sqrt(k * k - za2);

  double w2 = z * z / pow(g + n - fabs((double)k), 2);
  double d = 1. + a2 * w2;

  return -w2 / (2 * d) - (0.5 * a2 * w2 + 1. - sqrt(1. + a2 * w2)) * (c2 / d);
}

} // namespace AtomInfo

//******************************************************************************
struct NonRelSEConfig {
  int n;
  int l;
  int num;
  NonRelSEConfig(int in_n, int in_l, int in_num)
      : n(in_n), l(in_l), num(in_num) {}
  NonRelSEConfig() {}
  std::string symbol() {
    return std::to_string(n) + AtomInfo::l_symbol(l) + std::to_string(num);
  }

  bool ok() {
    if (l + 1 > n || l < 0 || n < 1 || num < 0 || num > 4 * l + 2)
      return false;
    return true;
  }

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

struct DiracSEnken { // name OK? too short?
  int n;
  int k;
  double en;
  DiracSEnken(int in_n, int in_k, double in_en = 0)
      : n(in_n), k(in_k), en(in_en){};
  DiracSEnken(){}; // never used? make above const!
};
