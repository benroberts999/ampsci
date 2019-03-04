#pragma once
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace ATI {

// Default values for A for each atom.
// Note: array index matches Z, so first entry is blank.
// Goes up to E120 (Z=120)
static const int A[121] = {
    0,   1,   4,   7,   9,   11,  12,  14,  16,  19,  20,  23,  24,  27,
    28,  31,  32,  35,  40,  39,  40,  45,  48,  51,  52,  55,  56,  59,
    59,  64,  65,  70,  73,  75,  79,  80,  84,  85,  88,  89,  91,  93,
    96,  97,  101, 103, 106, 108, 112, 115, 119, 122, 128, 127, 131, 133,
    137, 139, 140, 141, 144, 145, 150, 152, 157, 159, 162, 165, 167, 169,
    173, 175, 178, 181, 184, 186, 190, 192, 195, 197, 201, 204, 207, 209,
    209, 210, 222, 223, 226, 227, 232, 231, 238, 237, 244, 243, 247, 247,
    251, 252, 257, 258, 259, 262, 267, 270, 269, 270, 270, 278, 281, 281,
    285, 286, 289, 289, 293, 293, 294, 315, 320};

static const std::string atom_name_z[121] = {
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

constexpr int defaultA(int Z) { return (Z < 121 && Z > 0) ? A[Z] : 0; }

inline std::string atomicSymbol(int Z) {
  return (Z < 121 && Z > 0) ? atom_name_z[Z] : std::to_string(Z);
}

// Given an atomic symbol (H, He, etc.), will return Z
// Note: Symbol must be exact, including capitalisation
inline int get_z(const std::string at) {
  for (int z = 0; z < 121; z++)
    if (at == atom_name_z[z])
      return z;
  return std::stoi(at);
}

static const std::string spectroscopic_notation = "spdfghiklmnoqrtuv";
static const std::string Spectroscopic_Notation = "SPDFGHIKLMNOQRTUV";

// Short function that returns orbital term given l
inline std::string l_symbol(int l) {
  if (l < (int)spectroscopic_notation.length() && l >= 0)
    return spectroscopic_notation.substr(l, 1);
  else
    return "[" + std::to_string(l) + "]";
}

inline int symbol_to_l(std::string l_str) {
  // const char?
  for (int i = 0; i < (int)spectroscopic_notation.length(); i++)
    if (spectroscopic_notation.substr(i, 1) == l_str)
      return i;
  std::cerr << "\nFAIL ATI::69 Invalid l: " << l_str << "?\n";
  return -1;
}

constexpr int l_k(int ka) { return (abs(2 * ka + 1) - 1) / 2; }
constexpr int twoj_k(int ka) { return 2 * abs(ka) - 1; }
constexpr double j_k(int ka) { return 0.5 * twoj_k(ka); }

constexpr int indexFromKappa(int ka) {
  return (ka < 0) ? -2 * ka - 2 : 2 * ka - 1;
}

constexpr int kappaFromIndex(int i) {
  return (i % 2 == 0) ? -(i + 2) / 2 : (i + 1) / 2;
}

// const std::string core_config_list = {"He", "Ne", "Ar", "Kr", "Xe",
// "Rn", "Og", "Zn", "Cd", "Hg"};

// Shell configurations for Noble gasses (Group 8)
static const std::vector<int> core_He = {2};
static const std::vector<int> core_Ne = {2, 2, 6};
static const std::vector<int> core_Ar = {2, 2, 6, 2, 6};
static const std::vector<int> core_Kr = {2, 2, 6, 2, 6, 10, 2, 6};
static const std::vector<int> core_Xe = {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2, 6};
static const std::vector<int> core_Rn = {2,  2, 6, 2,  6, 10, 2, 6, 10,
                                         14, 2, 6, 10, 0, 0,  2, 6};
static const std::vector<int> core_Og = {
    2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 14, 0, 2, 6, 10, 0, 0, 0, 2, 6};

// Some other useful 'semi' full shells (transition)
static const std::vector<int> core_Zn = {2, 2, 6, 2, 6, 10, 2};
static const std::vector<int> core_Cd = {2, 2, 6, 2, 6, 10, 2, 6, 10, 0, 2};
static const std::vector<int> core_Hg = {2,  2,  6, 2, 6,  10, 2, 6,
                                         10, 14, 2, 6, 10, 0,  0, 2};

const std::vector<int> core_n = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5,
                                 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8,
                                 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9};
const std::vector<int> core_l = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4,
                                 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1,
                                 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8};

inline std::vector<int> getCoreConfig(const std::string ng) {
  if (ng == "He")
    return ATI::core_He;
  else if (ng == "Ne")
    return ATI::core_Ne;
  else if (ng == "Ar")
    return ATI::core_Ar;
  else if (ng == "Kr")
    return ATI::core_Kr;
  else if (ng == "Xe")
    return ATI::core_Xe;
  else if (ng == "Rn")
    return ATI::core_Rn;
  else if (ng == "Og")
    return ATI::core_Og;
  else if (ng == "Zn")
    return ATI::core_Zn;
  else if (ng == "Cd")
    return ATI::core_Cd;
  else if (ng == "Hg")
    return ATI::core_Hg;
  return {};
}

//******************************************************************************
inline double diracen(double z, double n, int k,
                      double alpha = 0.00729735256635) {
  double a2 = pow(alpha, 2);
  double c2 = 1. / pow(alpha, 2);
  double za2 = pow(alpha * z, 2);
  double g = sqrt(k * k - za2);

  double w2 = pow(z, 2) / pow(g + n - fabs((double)k), 2);
  double d = 1. + a2 * w2;

  return -w2 / (2 * d) - (a2 * w2 / 2 + 1. - sqrt(1. + a2 * w2)) * (c2 / d);
}

} // namespace ATI
