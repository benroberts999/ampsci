#pragma once
#include <array>
#include <string>
#include <utility>
#include <vector>

namespace AtomData {

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
  std::string name;
  Element(int inZ, const std::string &insymbol, int inA,
          const std::string &inname)
      : Z(inZ), symbol(insymbol), A(inA), name(inname) {}
};

// Default values for A for each atom.
// Goes up to E120 (Z=120)
static const std::vector<Element> periodic_table = {
    {1, "H", 1, "hydrogen"},
    {2, "He", 4, "helium"},
    {3, "Li", 7, "lithium"},
    {4, "Be", 9, "beryllium"},
    {5, "B", 11, "boron"},
    {6, "C", 12, "carbon"},
    {7, "N", 14, "nitrogen"},
    {8, "O", 16, "oxygen"},
    {9, "F", 19, "fluorine"},
    {10, "Ne", 20, "neon"},
    {11, "Na", 23, "sodium"},
    {12, "Mg", 24, "magnesium"},
    {13, "Al", 27, "aluminium"},
    {14, "Si", 28, "silicon"},
    {15, "P", 31, "phosphorus"},
    {16, "S", 32, "sulphur"},
    {17, "Cl", 35, "chlorine"},
    {18, "Ar", 40, "argon"},
    {19, "K", 39, "potassium"},
    {20, "Ca", 40, "calcium"},
    {21, "Sc", 45, "scandium"},
    {22, "Ti", 48, "titanium"},
    {23, "V", 51, "vanadium"},
    {24, "Cr", 52, "chromium"},
    {25, "Mn", 55, "manganese"},
    {26, "Fe", 56, "iron"},
    {27, "Co", 59, "cobalt"},
    {28, "Ni", 59, "nickel"},
    {29, "Cu", 64, "copper"},
    {30, "Zn", 65, "zinc"},
    {31, "Ga", 70, "gallium"},
    {32, "Ge", 73, "germanium"},
    {33, "As", 75, "arsenic"},
    {34, "Se", 79, "selenium"},
    {35, "Br", 80, "bromine"},
    {36, "Kr", 84, "krypton"},
    {37, "Rb", 85, "rubidium"},
    {38, "Sr", 88, "strontium"},
    {39, "Y", 89, "yttrium"},
    {40, "Zr", 91, "zirconium"},
    {41, "Nb", 93, "niobium"},
    {42, "Mo", 96, "molybdenum"},
    {43, "Tc", 97, "technetium"},
    {44, "Ru", 101, "ruthenium"},
    {45, "Rh", 103, "rhodium"},
    {46, "Pd", 106, "palladium"},
    {47, "Ag", 108, "silver"},
    {48, "Cd", 112, "cadmium"},
    {49, "In", 115, "indium"},
    {50, "Sn", 119, "tin"},
    {51, "Sb", 122, "antimony"},
    {52, "Te", 128, "tellurium"},
    {53, "I", 127, "iodine"},
    {54, "Xe", 131, "xenon"},
    {55, "Cs", 133, "cesium"},
    {56, "Ba", 137, "barium"},
    {57, "La", 139, "lanthanum"},
    {58, "Ce", 140, "cerium"},
    {59, "Pr", 141, "praseodymium"},
    {60, "Nd", 144, "neodymium"},
    {61, "Pm", 145, "promethium"},
    {62, "Sm", 150, "samarium"},
    {63, "Eu", 152, "europium"},
    {64, "Gd", 157, "gadolinium"},
    {65, "Tb", 159, "terbium"},
    {66, "Dy", 162, "dysprosium"},
    {67, "Ho", 165, "holmium"},
    {68, "Er", 167, "erbium"},
    {69, "Tm", 169, "thulium"},
    {70, "Yb", 173, "ytterbium"},
    {71, "Lu", 175, "lutetium"},
    {72, "Hf", 178, "hafnium"},
    {73, "Ta", 181, "tantalum"},
    {74, "W", 184, "tungsten"},
    {75, "Re", 186, "rhenium"},
    {76, "Os", 190, "osmium"},
    {77, "Ir", 192, "iridium"},
    {78, "Pt", 195, "platinum"},
    {79, "Au", 197, "gold"},
    {80, "Hg", 201, "mercury"},
    {81, "Tl", 204, "thallium"},
    {82, "Pb", 207, "lead"},
    {83, "Bi", 209, "bismuth"},
    {84, "Po", 209, "polonium"},
    {85, "At", 210, "astatine"},
    {86, "Rn", 222, "radon"},
    {87, "Fr", 223, "francium"},
    {88, "Ra", 226, "radium"},
    {89, "Ac", 227, "actinium"},
    {90, "Th", 232, "thorium"},
    {91, "Pa", 231, "protactinium"},
    {92, "U", 238, "uranium"},
    {93, "Np", 237, "neptunium"},
    {94, "Pu", 244, "plutonium"},
    {95, "Am", 243, "americium"},
    {96, "Cm", 247, "curium"},
    {97, "Bk", 247, "berkelium"},
    {98, "Cf", 251, "californium"},
    {99, "Es", 252, "einsteinium"},
    {100, "Fm", 257, "fermium"},
    {101, "Md", 258, "mendelevium"},
    {102, "No", 259, "nobelium"},
    {103, "Lr", 262, "lawrencium"},
    {104, "Rf", 267, "rutherfordium"},
    {105, "Db", 270, "dubnium"},
    {106, "Sg", 269, "seaborgium"},
    {107, "Bh", 270, "bohrium"},
    {108, "Hs", 270, "hassium"},
    {109, "Mt", 278, "meitnerium"},
    {110, "Ds", 281, "darmstadtium"},
    {111, "Rg", 281, "roentgenium"},
    {112, "Cn", 285, "copernicium"},
    {113, "Nh", 286, "nihonium"},
    {114, "Fl", 289, "flerovium"},
    {115, "Mc", 289, "moscovium"},
    {116, "Lv", 293, "livermorium"},
    {117, "Ts", 293, "tennessine"},
    {118, "Og", 294, "oganesson"},
    {119, "E119", 315, "ununennium"},
    {120, "E120", 320, "unbinilium"}};

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

} // namespace AtomData
