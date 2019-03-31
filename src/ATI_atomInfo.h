#pragma once
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace ATI {

// Default values for A for each atom.
// Note: array index matches Z, so first entry is blank.
// Goes up to E120 (Z=120)
// static const int
// XXX Make this a pair? Or a struct? (then, can put other default nuclear
// defaults inside as well!)
static const std::size_t MAX_Z = 121;
static const std::array<int, MAX_Z> A = {
    0,   1,   4,   7,   9,   11,  12,  14,  16,  19,  20,  23,  24,  27,
    28,  31,  32,  35,  40,  39,  40,  45,  48,  51,  52,  55,  56,  59,
    59,  64,  65,  70,  73,  75,  79,  80,  84,  85,  88,  89,  91,  93,
    96,  97,  101, 103, 106, 108, 112, 115, 119, 122, 128, 127, 131, 133,
    137, 139, 140, 141, 144, 145, 150, 152, 157, 159, 162, 165, 167, 169,
    173, 175, 178, 181, 184, 186, 190, 192, 195, 197, 201, 204, 207, 209,
    209, 210, 222, 223, 226, 227, 232, 231, 238, 237, 244, 243, 247, 247,
    251, 252, 257, 258, 259, 262, 267, 270, 269, 270, 270, 278, 281, 281,
    285, 286, 289, 289, 293, 293, 294, 315, 320};

static const std::array<std::string, MAX_Z> atom_name_z = {
    "0",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",    "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",   "Sc",
    "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",   "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo",   "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",    "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",   "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re",   "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",   "Fr",
    "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk",   "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",   "Mt",
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "E119", "E120"};

constexpr int defaultA(int Z) { return (Z < (int)MAX_Z && Z > 0) ? A[Z] : 0; }

inline std::string atomicSymbol(int Z) {
  return (Z < (int)MAX_Z && Z > 0) ? atom_name_z[Z] : std::to_string(Z);
}

// Given an atomic symbol (H, He, etc.), will return Z
// Note: Symbol must be exact, including capitalisation
inline int get_z(const std::string &at) {
  for (int z = 0; z < (int)MAX_Z; z++) {
    if (at == atom_name_z[z])
      return z;
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
    std::cerr << "\nFAIL ATI::69 Invalid l: " << l_str << "?\n";
  }
  return l;
}

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

// XXX This is a bad/limiting solution:
// two issues: a) retrive (n,l) from index
// b) parsing the input string [not big deal, n typically <7]
const std::array<int, 45> core_n = {
    1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7,
    7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9};
const std::array<int, 45> core_l = {
    0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1,
    2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8};

// Note: this requires that all Nobel Gasses are listed FIRST, in order
// (Assumed by "niceCoreOutput" function that this matches nobelGasses
static const std::array<std::pair<std::string, std::string>, 11> nobelGasses = {
    std::make_pair("[He]", "1s2"), /**/ //
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
                           "4f14,5d10,6s2,6p6,5f14,6d10,7s2")};

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

} // namespace ATI
