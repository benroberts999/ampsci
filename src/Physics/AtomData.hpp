#pragma once
#include <array>
#include <string>
#include <vector>

//******************************************************************************
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
    return (num < filling) ? double(num) / double(filling) : 1;
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

//******************************************************************************
struct DiracSEnken { // name OK? too short?
  int n;
  int k;
  double en;
  DiracSEnken(int in_n = 0, int in_k = 0, double in_en = 0)
      : n(in_n), k(in_k), en(in_en){};
};

//******************************************************************************
//******************************************************************************
namespace AtomData {

int defaultA(int Z);

std::string atomicSymbol(int Z);
std::string atomicName(int Z);

inline int get_z(int z) { return z; }
int get_z(const std::string &at);

std::string l_symbol(int l);

int symbol_to_l(const std::string &l_str);

std::string coreConfig(const std::string &in_ng);

std::string niceCoreOutput(const std::string &full_core);

double diracen(double z, double n, int k, double alpha = 0.00729735256635);

std::vector<NonRelSEConfig> core_parser(const std::string &str_core_in);

std::string guessCoreConfigStr(const int total_core_electrons);
std::vector<NonRelSEConfig> core_guess(const int total_core_electrons);

std::vector<DiracSEnken> listOfStates_nk(const std::string &in_list);

void printTable();

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

} // namespace AtomData
