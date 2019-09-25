#pragma once
#include <array>
#include <string>
#include <utility>
#include <vector>

namespace AtomInfo {

static const std::string spectroscopic_notation = "spdfghiklmnoqrtuvwxyzabc";
static const std::string Spectroscopic_Notation = "SPDFGHIKLMNOQRTUVWXYZABC";

const std::string filling_order =
    "1s2,2s2,2p6,3s2,3p6,4s2,3d10,4p6,5s2,4d10,5p6,6s2,4f14,5d10,"
    "6p6,7s2,5f14,6d10,7p6,8s2,6f14,7d10,8p6,9s2";

//******************************************************************************
struct Element {
  int Z;
  std::string symbol;
  int A;
  Element(int inZ, const std::string &insymbol, int inA)
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

// Note: this requires that all Nobel Gasses are listed FIRST, in order
// (Assumed by "niceCoreOutput" function that this matches nobelGasses)
// Nothing functional will break if not, just won't necisarily find 'nicest'
// output format
static const std::array<std::pair<std::string, std::string>, 13> nobelGasses = {
    std::make_pair("[]", ""),      //
    std::make_pair("[He]", "1s2"), //
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

} // namespace AtomInfo
