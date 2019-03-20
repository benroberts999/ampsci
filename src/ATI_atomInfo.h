#pragma once
#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ATI {

// Default values for A for each atom.
// Note: array index matches Z, so first entry is blank.
// Goes up to E120 (Z=120)
// static const int
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

// static const std::string atom_name_z[121] =
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
inline int get_z(const std::string at) {
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

static const std::string spectroscopic_notation = "spdfghiklmnoqrtuvwxyz";
static const std::string Spectroscopic_Notation = "SPDFGHIKLMNOQRTUVWXYZ";

// Short function that returns orbital term given l
inline std::string l_symbol(int l) {
  if (l < (int)spectroscopic_notation.length() && l >= 0)
    return spectroscopic_notation.substr(l, 1);
  else
    return "[" + std::to_string(l) + "]";
}

inline int symbol_to_l(std::string l_str) {
  // const char?
  for (int i = 0; i < (int)spectroscopic_notation.length(); i++) {
    if (spectroscopic_notation.substr(i, 1) == l_str)
      return i;
  }
  int l = -1;
  try {
    // Can work if given an int as a string:
    l = std::stoi(l_str);
  } catch (...) {
    std::cerr << "\nFAIL ATI::69 Invalid l: " << l_str << "?\n";
  }
  return l;
}

inline int l_k(int ka) { return (ka > 0) ? ka : -ka - 1; }
inline int twoj_k(int ka) { return 2 * abs(ka) - 1; }
inline double j_k(int ka) { return abs(ka) - 0.5; }
inline int parity_k(int ka) {
  return (ka % 2 == 0) ? ((ka > 0) ? 1 : -1) : ((ka < 0) ? 1 : -1);
}

constexpr int indexFromKappa(int ka) {
  return (ka < 0) ? -2 * ka - 2 : 2 * ka - 1;
}

constexpr int kappaFromIndex(int i) {
  return (i % 2 == 0) ? -(i + 2) / 2 : (i + 1) / 2;
}

// XXX This is a bad/limiting solution:
const std::vector<int> core_n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5,
                                 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8,
                                 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9};
const std::vector<int> core_l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4,
                                 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1,
                                 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8};

inline std::string coreConfig(const std::string &ng) {
  // NOTE: must return SAME string if none found (or, replace?)
  if (ng == "He")
    return "1s2";
  else if (ng == "Ne")
    return coreConfig("He") + ",2s2,2p6";
  else if (ng == "Ar")
    return coreConfig("Ne") + ",3s2,3p6";
  else if (ng == "Kr")
    return coreConfig("Ar") + ",3d10,4s2,4p6";
  else if (ng == "Xe")
    return coreConfig("Kr") + ",4d10,5s2,5p6";
  else if (ng == "Rn")
    return coreConfig("Xe") + ",4f14,5d10,6s2,6p6";
  else if (ng == "Og")
    return coreConfig("Rn") + ",5f14,6d10,7s2,7p6";
  else if (ng == "Zn")
    return coreConfig("Ar") + ",3d10,4s2";
  else if (ng == "Cd")
    return coreConfig("Kr") + ",4d10,5s2";
  else if (ng == "Hg")
    return coreConfig("Xe") + ",4f14,5d10,6s2";
  else if (ng == "Cn")
    return coreConfig("Rn") + ",5f14,6d10,7s2";
  else if (ng == "Yb")
    return coreConfig("Xe") + ",4f14,6s2";
  else if (ng == "No")
    return coreConfig("Rn") + ",5f14,7s2";
  else
    return ng;
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
