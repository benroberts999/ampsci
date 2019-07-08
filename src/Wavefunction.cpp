// class Wavefunction::
#include "Wavefunction.hpp"
#include "ADAMS_bound.hpp"
#include "AtomInfo.hpp"
#include "DiracSpinor.hpp"
#include "Grid.hpp"
#include "Nuclear.hpp"
#include "PhysConst_constants.hpp"
#include <algorithm> //for sort
#include <cmath>
#include <locale> //isdigit
#include <sstream>
#include <string>
#include <vector>

//******************************************************************************
Wavefunction::Wavefunction(int in_z, int in_a, int in_ngp, double rmin,
                           double rmax, double var_alpha)
    : rgrid(rmin, rmax, (std::size_t)in_ngp, GridType::loglinear, 3.5),
      m_alpha(PhysConst::alpha * var_alpha), m_Z(in_z),
      m_A((in_a < 0) ? AtomInfo::defaultA(m_Z) : in_a) {
  // Make Vnuc const?
  if (m_A > 15) {
    formNuclearPotential(Nuclear::Type::Fermi);
  } else if (m_A > 0) {
    formNuclearPotential(Nuclear::Type::spherical);
  } else {
    formNuclearPotential(Nuclear::Type::zero);
  }

  if (m_Z * m_alpha > 1) {
    std::cerr << "Alpha too large, Za=" << m_Z * m_alpha << "\n";
    std::abort();
  }
}

//******************************************************************************
void Wavefunction::solveDirac(DiracSpinor &psi, double e_a,
                              const std::vector<double> &vex,
                              int log_dele_or) const
// Uses ADAMS::solveDBS to solve Dirac Eqn for local potential (Vnuc + Vdir)
// If no e_a is given, will use the existing one!
// (Usually, a better guess should be given, using P.T.)
// Note: optionally takes in exchange potential! (see overloaded above)
// Note: Uses the "dodgy" re-scaled exchange potenital:
// Vex\psi_a = sum_b vex_a psi_b -> [sum_b vex_a (psi_b/psi_a)] psi_a
// so, here, vex = [sum_b vex_a (psi_b/psi_a)]
// This is not ideal..
{
  std::vector<double> v_a = vnuc;
  if (vdir.size() != 0) {
    for (std::size_t i = 0; i < rgrid.ngp; i++) {
      v_a[i] += vdir[i];
    }
  }
  if (vex.size() != 0) {
    for (std::size_t i = 0; i < rgrid.ngp; i++) {
      v_a[i] += vex[i];
    }
  }
  if (e_a != 0) {
    psi.en = e_a;
  } else if (psi.en == 0) {
    psi.en = enGuessVal(psi.n, psi.k);
  }
  ADAMS::solveDBS(psi, v_a, rgrid, m_alpha, log_dele_or);
}

//------------------------------------------------------------------------------
void Wavefunction::solveDirac(DiracSpinor &psi, double e_a,
                              int log_dele_or) const
// Overloaded version; see above
// This one doesn't have exchange potential
{
  std::vector<double> empty_vec;
  return solveDirac(psi, e_a, empty_vec, log_dele_or);
}

// //******************************************************************************
// static std::vector<NonRelSEConfig> core_parser(const std::string
// &str_core_in)
// // Heler function for below.
// // Move to Atom Info ?
// {
//   // If there's a 'Noble-Gas' term, replace it with full config
//   // Otherwise, 'first-term' remains unchanges
//   auto found = str_core_in.find(",");
//   if (found > str_core_in.length())
//     found = str_core_in.length();
//   auto first_term = str_core_in.substr(0, found);
//   auto rest = str_core_in.substr(found);
//   auto str_core = AtomInfo::coreConfig(first_term) + rest;
//
//   // Move comma-seperated string into an array (vector)
//   std::vector<std::string> term_str_list;
//   {
//     std::stringstream ss(str_core);
//     while (ss.good()) {
//       std::string substr;
//       getline(ss, substr, ',');
//       term_str_list.push_back(substr);
//     }
//   }
//
//   std::vector<NonRelSEConfig> core_configs;
//   for (const auto &term : term_str_list) {
//     bool term_ok = true;
//     std::size_t l_position = 0;
//     for (const auto &c : term) {
//       if (!std::isdigit(c))
//         break;
//       ++l_position;
//     }
//     int n{0}, num{0}, l{-1};
//     try {
//       n = std::stoi(term.substr(0, l_position - 0));
//       num = std::stoi(term.substr(l_position + 1));
//       if (l_position == term.size())
//         throw;
//       l = AtomInfo::symbol_to_l(term.substr(l_position, 1));
//     } catch (...) {
//       term_ok = false;
//     }
//     NonRelSEConfig new_config(n, l, num);
//
//     if (!term_ok || n <= 0) {
//       std::cout << "Problem with core: " << str_core_in << "\n";
//       std::cerr << "invalid core term: " << term << "\n";
//       std::abort();
//     }
//
//     if (num == 0)
//       continue;
//     auto ia = std::find(core_configs.begin(), core_configs.end(),
//     new_config); if (ia == core_configs.end()) {
//       core_configs.push_back(new_config);
//     } else {
//       *ia += new_config;
//     }
//   }
//   return core_configs;
// }

//******************************************************************************
void Wavefunction::determineCore(const std::string &str_core_in)
// Takes in a string list for the core configuration, outputs an int list
// Takes in previous closed shell (noble), + 'rest' (or just the rest)
// E.g:
//   Core of Cs: Xe
//   Core of Gold: Xe 4f14 5d10
// 'rest' is in form nLm : n=n, L=l, m=number of electrons in that nl shell.
{

  m_core_configs = AtomInfo::core_parser(str_core_in);

  bool bad_core = false;
  m_core_string = "";
  for (auto config : m_core_configs) {
    if (!config.ok()) {
      m_core_string += " **";
      bad_core = true;
    }
    m_core_string += config.symbol();
    if (config != m_core_configs.back())
      m_core_string += ",";
  }

  if (bad_core) {
    std::cout << "Problem with core: " << str_core_in << " = \n";
    std::cout << m_core_string << "\n";
    std::cout << "In this house, we obey the Pauli exclusion principle!\n";
    std::abort();
  }

  // Count number of electrons in the core
  num_core_electrons = 0;
  for (const auto &config : m_core_configs) {
    num_core_electrons += config.num;
  }

  if (num_core_electrons > m_Z) {
    std::cout << "Problem with core: " << str_core_in << "\n";
    std::cout << "= " << m_core_string << "\n";
    std::cout << "= " << AtomInfo::niceCoreOutput(m_core_string) << "\n";
    std::cout << "Too many electrons: N_core=" << num_core_electrons
              << ", Z=" << m_Z << "\n";
    std::abort();
  }

  return;
}

//******************************************************************************
bool Wavefunction::isInCore(int n, int k) const
// Checks if given state is in the core.
{
  for (auto &phi : core_orbitals) {
    if (n == phi.n && k == phi.k)
      return true;
  }
  return false;
}

//******************************************************************************
std::size_t Wavefunction::getStateIndex(int n, int k, bool &is_valence) const {

  is_valence = false;
  for (const auto p_orbs : {&core_orbitals, &valence_orbitals}) {
    std::size_t count = 0;
    for (const auto &phi : *p_orbs) {
      if (n == phi.n && k == phi.k)
        return count;
      ++count;
    }
    is_valence = true;
  }
  std::cerr << "\nFAIL 290 in WF: Couldn't find state n,k=" << n << "," << k
            << "\n";
  std::abort();
}
std::size_t Wavefunction::getStateIndex(const DiracSpinor &psi,
                                        bool &is_valence) const {
  return getStateIndex(psi.n, psi.k, is_valence);
}

//******************************************************************************
int Wavefunction::maxCore_n(int ka) const
// Returns the largest n for states with kappa = ka in the core
// Note: ka is optional input; if none given, will be 0 (& not used)
// (used for energy guesses)
// Note: if you give it l instead of kappa, still works!
{
  int max_n = 0;
  for (auto &phi : core_orbitals) {
    if (phi.k != ka && ka != 0)
      continue;
    if (phi.n > max_n)
      max_n = phi.n;
  }
  return max_n;
}

//******************************************************************************
void Wavefunction::solveInitialCore(std::string str_core, int log_dele_or)
// Solves the Dirac eqn for each state in the core
// Only for local potential (direct part)
// HartreeFockClass.cpp has routines for Hartree Fock
{
  if (core_orbitals.size() > 0) {
    std::cerr << "Fail 254 in Wavefunction:solveInitialCore: States "
                 "already exist! "
              << core_orbitals.size() << "\n";
    std::abort();
  }

  determineCore(str_core);

  for (const auto &config : m_core_configs) {
    int num = config.num;
    if (num == 0)
      continue;
    int n = config.n;
    int l = config.l;
    double en_a = enGuessCore(n, l);
    int k1 = l; // j = l-1/2
    if (k1 != 0) {
      core_orbitals.emplace_back(DiracSpinor{n, k1, rgrid});
      solveDirac(core_orbitals.back(), en_a, log_dele_or);
      core_orbitals.back().occ_frac = double(num) / (4 * l + 2);
      en_a = 0.95 * core_orbitals.back().en;
      if (en_a > 0)
        en_a = enGuessCore(n, l);
    }
    int k2 = -(l + 1); // j=l+1/2
    core_orbitals.emplace_back(DiracSpinor{n, k2, rgrid});
    solveDirac(core_orbitals.back(), en_a, log_dele_or);
    core_orbitals.back().occ_frac = double(num) / (4 * l + 2);
  }
}

//******************************************************************************
void Wavefunction::solveNewValence(int n, int k, double en_a, int log_dele_or)
// Update to take a list ok nken's ?
{
  valence_orbitals.emplace_back(DiracSpinor{n, k, rgrid});

  // Solve local dirac Eq:
  auto &psi = valence_orbitals.back();
  if (en_a == 0)
    en_a = enGuessVal(n, k);
  solveDirac(psi, en_a, log_dele_or);
}

//******************************************************************************
void Wavefunction::orthonormaliseOrbitals(std::vector<DiracSpinor> &in_orbs,
                                          int num_its)
// Note: this function is static
// Forces ALL orbitals to be orthogonal to each other, and normal
// Note: workes best if run twice!
// |a> ->  |a> - \sum_{b!=a} |b><b|a>
// Then:
// |a> -> |a> / <a|a>
// c_ba = c_ab = <a|b>
// num_its is optional parameter. Repeats that many times!
// Note: I force all orthog to each other - i.e. double count.
// {force <2|1>=0 and then <1|2>=0}
// Would be 2x faster not to do this
//  - but that would treat some orbitals special!
// Hence factor of 0.5
// Note: For HF, should never be called after core is frozen!
//
// Note: This allows wfs to extend past pinf!
// ==> This causes the possible orthog issues..
{
  auto Ns = in_orbs.size();
  auto ngp = in_orbs.front().p_rgrid->ngp;

  std::vector<std::vector<double>> c_ab(Ns, std::vector<double>(Ns));
  // Calculate c_ab = <a|b>  [only for b>a -- symmetric]
  for (std::size_t a = 0; a < Ns; a++) {
    for (auto b = a + 1; b < Ns; b++) {
      if (in_orbs[a].k != in_orbs[b].k)
        continue;
      c_ab[a][b] = 0.5 * (in_orbs[a] * in_orbs[b]);
    }
  }

  // Orthogonalise orbitals:
  for (std::size_t a = 0; a < Ns; a++) {
    for (std::size_t b = 0; b < Ns; b++) {
      if (in_orbs[a].k != in_orbs[b].k)
        continue;
      if (a == b)
        continue;
      // c_ab = c_ba : only calc'd half:
      double cab = (a < b) ? c_ab[a][b] : c_ab[b][a];
      for (std::size_t ir = 0; ir < ngp; ir++) {
        in_orbs[a].f[ir] -= cab * in_orbs[b].f[ir];
        in_orbs[a].g[ir] -= cab * in_orbs[b].g[ir];
      }
    }
  }

  for (auto &psi : in_orbs) {
    psi.normalise();
  }

  // If necisary: repeat
  if (num_its > 1)
    orthonormaliseOrbitals(in_orbs, num_its - 1);
}

//******************************************************************************
void Wavefunction::orthonormaliseWrtCore(DiracSpinor &psi_v) const
// Force given orbital to be orthogonal to all core orbitals
// [After the core is 'frozen', don't touch core orbitals!]
// |v> --> |v> - sum_c |c><c|v>
// note: here, c denotes core orbitals
{
  // Orthogonalise:
  for (const auto &psi_c : core_orbitals) {
    if (psi_v.k != psi_c.k)
      continue;
    const double Avc = psi_v * psi_c;
    for (std::size_t ir = 0; ir < psi_v.pinf; ir++) {
      psi_v.f[ir] -= Avc * psi_c.f[ir];
      psi_v.g[ir] -= Avc * psi_c.g[ir];
    }
  }
  psi_v.normalise();
}

//******************************************************************************
double Wavefunction::enGuessCore(int n, int l) const
// Private
// Energy guess for core states. Not perfect, good enough
// num_el_below = total electrons BELOW
// num = num electrons in THIS shell
{
  int num_el_below = 0;
  int num_el_this = 0;
  for (const auto &config : m_core_configs) {
    if (l == config.l && n == config.l) {
      num_el_this = config.num;
      break;
    }
    num_el_below += config.num;
  }

  // effective Z (for energy guess) -- not perfect!
  double Zeff = (m_Z - num_el_below - num_el_this);
  if (l == 1) {
    Zeff = 1. + (m_Z - num_el_below - 0.5 * num_el_this);
  } else if (l == 2) {
    Zeff = 1. + (m_Z - num_el_below - 0.5 * num_el_this);
  }
  if (Zeff < 1.0) {
    Zeff = 1.;
  }

  double en_a = -0.5 * pow(Zeff / n, 2);
  if (n > 1) {
    en_a *= 0.5;
  }

  if (n == maxCore_n() - 1) {
    en_a *= 1.25;
  } else if (n == maxCore_n()) {
    en_a *= 1.5;
  } else if (n == maxCore_n() + 1 && l == 0 && num_el_this == 2) {
    en_a = -0.5 * pow(Zeff, 2);
  }

  return en_a;
}

//******************************************************************************
double Wavefunction::enGuessVal(int n, int ka) const
// Energy guess for valence states. Not perfect, good enough
{
  int maxn = maxCore_n();
  int l = AtomInfo::l_k(ka);
  int dn = n - maxn;
  double neff = 1. + dn;
  double x = 1;
  if (maxn < 4)
    x = 0.25;
  if (l == 1)
    neff += 0.5 * x;
  if (l == 2)
    neff += 2. * pow(x, 0.5);
  if (l >= 3)
    neff += 4. * x;
  return -0.5 / pow(neff, 2);
}

//******************************************************************************
void Wavefunction::formNuclearPotential(Nuclear::Type nucleus_type, double rc,
                                        double t) {
  vnuc.clear();
  switch (nucleus_type) {
  case Nuclear::Type::Fermi:
    if (t == 0)
      t = Nuclear::approximate_t_skin(m_A);
    if (rc == 0) {
      auto rrms = Nuclear::find_rrms(m_Z, m_A);
      if (rrms == 0)
        rrms = Nuclear::approximate_r_rms(m_A);
      rc = Nuclear::c_hdr_formula_rrms_t(rrms, t);
    }
    vnuc = Nuclear::fermiNuclearPotential(m_Z, t, rc, rgrid.r);
    break;
  case Nuclear::Type::spherical:
    if (rc == 0) {
      // note: still called rc, but is r_N here!
      rc = Nuclear::find_rrms(m_Z, m_A);
      if (rc == 0)
        rc = Nuclear::approximate_r_rms(m_A);
    }
    t = 0;
    vnuc = Nuclear::sphericalNuclearPotential(m_Z, rc, rgrid.r);
    break;
  case Nuclear::Type::zero:
    rc = 0;
    t = 0;
    vnuc = Nuclear::sphericalNuclearPotential(m_Z, 0., rgrid.r);
    break;
  default:
    std::cerr << "\nFail WF:755 - invalid nucleus type?\n";
  }
  m_c = rc;
  m_t = t;
}

//******************************************************************************
std::string Wavefunction::nuclearParams() const {
  std::ostringstream output;
  if (m_c == 0 && m_t == 0) {
    output << "Zero-size nucleus";
  } else if (m_t == 0) {
    output << "Spherical nucleus; r_rms = " << m_c;
  } else {
    output << "Fermi nucleus; r_rms = " << Nuclear::rrms_formula_c_t(m_c, m_t)
           << ", c=" << m_c << ", t=" << m_t;
  }
  return output.str();
}

//******************************************************************************
std::vector<std::size_t>
Wavefunction::sortedEnergyList(const std::vector<DiracSpinor> &tmp_orbs,
                               bool do_sort) const
// Outouts a list of integers corresponding to the states
// sorted by energy (lowest energy first)
{

  using DoubleInt = std::pair<double, std::size_t>;
  std::vector<DoubleInt> t_en;
  for (std::size_t i = 0; i < tmp_orbs.size(); i++) {
    t_en.emplace_back(tmp_orbs[i].en, i);
  }

  // Sort list of Pairs by first element in the pair:
  auto compareThePair = [](const DoubleInt &di1, const DoubleInt &di2) {
    return di1.first < di2.first;
  };

  if (do_sort)
    std::sort(t_en.begin(), t_en.end(), compareThePair);

  // overwrite list with sorted list
  std::vector<std::size_t> sorted_list;
  std::for_each(t_en.begin(), t_en.end(),
                [&](const DoubleInt &el) { sorted_list.push_back(el.second); });

  return sorted_list;
}

//******************************************************************************
void Wavefunction::printCore(bool sorted) const
// prints core orbitals
{
  int Zion = Znuc() - Ncore();
  std::cout << "Core: " << coreConfiguration_nice() << " (V^N";
  if (Zion != 0)
    std::cout << "-" << Zion;
  std::cout << ")\n";
  std::cout << "     state   k   Rinf its   eps       En (au)      En (/cm)\n";

  auto index_list = sortedEnergyList(core_orbitals, sorted);
  for (auto i : index_list) {
    auto &phi = core_orbitals[i];
    double r_inf = rinf(phi);
    printf("%2i) %7s %2i  %5.1f %2i  %5.0e %13.7f %13.1f", int(i),
           phi.symbol().c_str(), phi.k, r_inf, phi.its, phi.eps, phi.en,
           phi.en *PhysConst::Hartree_invcm);
    if (phi.occ_frac < 1.0) {
      printf("     (%4.2f)\n", phi.occ_frac);
    } else {
      std::cout << "\n";
    }
  }
}

//******************************************************************************
void Wavefunction::printValence(
    bool sorted, const std::vector<DiracSpinor> &in_orbitals) const
// prints valence orbitals
{
  auto tmp_orbs = (in_orbitals.size() == 0) ? valence_orbitals : in_orbitals;

  if (tmp_orbs.size() == 0)
    return;

  std::cout << "Val: state   "
            << "k   Rinf its   eps       En (au)      En (/cm)   En (/cm)\n";

  // Find lowest valence energy:
  double e0 = 0;
  for (auto &phi : tmp_orbs) {
    if (phi.en < e0)
      e0 = phi.en;
  }

  auto index_list = sortedEnergyList(tmp_orbs, sorted);
  for (auto i : index_list) {
    auto &phi = tmp_orbs[i];
    double r_inf = rinf(phi);
    printf("%2i) %7s %2i  %5.1f %2i  %5.0e %13.7f %13.1f", int(i),
           phi.symbol().c_str(), phi.k, r_inf, phi.its, phi.eps, phi.en,
           phi.en *PhysConst::Hartree_invcm);
    printf(" %10.2f\n", (phi.en - e0) * PhysConst::Hartree_invcm);
  }
}

//******************************************************************************
std::vector<DiracSEnken>
Wavefunction::listOfStates_nk(int num_val, int la, int lb, bool skip_core) const
// Creates a list of states (usually valence states to solve for)
// In form {{n,ka},...}
// Outputs n and l, from min-max l (or from 0-> max l)
// Calulates num_val different n states (for each l)
// This will be number of states above the core (skips states which are in the
// core)
{
  std::vector<DiracSEnken> lst;
  auto l_min = la;
  auto l_max = lb;
  if (lb == 0) {
    l_min = 0;
    l_max = la;
  }

  auto min_ik = AtomInfo::indexFromKappa(-l_min - 1);
  auto max_ik = AtomInfo::indexFromKappa(-l_max - 1);
  for (int ik = min_ik; ik <= max_ik; ik++) {
    auto k = AtomInfo::kappaFromIndex(ik);
    auto l = AtomInfo::l_k(k);
    auto n_min = l + 1;
    for (int n = n_min, count = 0; count < num_val; n++) {
      if (isInCore(n, k) && skip_core)
        continue;
      lst.emplace_back(n, k);
      ++count;
    }
  }
  return lst;
}
