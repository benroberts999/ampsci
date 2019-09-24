#include "Dirac/Wavefunction.hpp"
#include "Adams/Adams_bound.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomInfo.hpp"
#include "Physics/Nuclear.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <algorithm> //for sort
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

//******************************************************************************
Wavefunction::Wavefunction(int z, const GridParameters &gridparams,
                           const Nuclear::Parameters &nuc_params,
                           double var_alpha)
    : rgrid({gridparams}),                                        //
      m_alpha(PhysConst::alpha * var_alpha),                      //
      m_Z(z), m_A(nuc_params.a),                                  //
      m_nuc_params(nuc_params),                                   //
      vnuc(Nuclear::formPotential(nuc_params, m_Z, m_A, rgrid.r)) //
//
{
  if (m_Z * m_alpha > 1) {
    std::cerr << "Alpha too large: Z*alpha=" << m_Z * m_alpha << "\n";
    std::abort();
  }
}

//******************************************************************************
void Wavefunction::solveDirac(DiracSpinor &psi, double e_a,
                              const std::vector<double> &vex,
                              int log_dele_or) const
// Uses Adams::solveDBS to solve Dirac Eqn for local potential (Vnuc + Vdir)
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
  Adams::solveDBS(psi, v_a, rgrid, m_alpha, log_dele_or);
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

//******************************************************************************
void Wavefunction::determineCore(std::string str_core_in)
// Takes in a string list for the core configuration, outputs an int list
// Takes in previous closed shell (noble), + 'rest' (or just the rest)
// E.g:
//   Core of Cs: Xe
//   Core of Gold: Xe 4f14 5d10
// 'rest' is in form nLm : n=n, L=l, m=number of electrons in that nl shell.
{

  // Check if integer; if so, V^N-M, '-M' is input integer.
  // Use 'guess' for core
  auto first_char = str_core_in.substr(0, 1);
  if (first_char == "0" || first_char == "-") {
    try {
      auto m = std::stoi(str_core_in);
      str_core_in = AtomInfo::guessCoreConfigStr(m_Z + m);
    } catch (...) {
    }
  }

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
void Wavefunction::hartreeFockCore(HFMethod method, const std::string &in_core,
                                   double eps_HF, double h_d, double g_t) {
  // XXX Update this (and HF) so that it doesn't re-Create m_pHF
  // AND, so that can re-run core!
  m_pHF = std::make_unique<HartreeFock>(
      HartreeFock(method, *this, in_core, eps_HF, h_d, g_t));
}

//******************************************************************************
auto Wavefunction::coreEnergyHF() const {
  if (!m_pHF) {
    std::cerr << "WARNING 62: Cant call coreEnergyHF before hartreeFockCore\n";
    return 0.0;
  }
  return m_pHF->calculateCoreEnergy();
}

//******************************************************************************
void Wavefunction::hartreeFockValence(const std::string &in_valence_str) {
  if (!m_pHF) {
    std::cerr << "WARNING 62: Cant call hartreeFockValence before "
                 "hartreeFockCore\n";
    return;
  }
  auto val_lst = AtomInfo::listOfStates_nk(in_valence_str);
  for (const auto &nk : val_lst) {
    if (!isInCore(nk.n, nk.k))
      m_pHF->solveNewValence(nk.n, nk.k);
  }
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
bool Wavefunction::isInCore(const DiracSpinor &phi) const {
  return isInCore(phi.n, phi.k);
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
  // std::abort();
  return std::max(core_orbitals.size(), valence_orbitals.size()); // invalid
}
std::size_t Wavefunction::getStateIndex(const DiracSpinor &psi,
                                        bool &is_valence) const {
  return getStateIndex(psi.n, psi.k, is_valence);
}
//******************************************************************************
const DiracSpinor &Wavefunction::getState(int n, int k,
                                          bool &is_valence) const {
  auto index = getStateIndex(n, k, is_valence);
  return is_valence ? valence_orbitals[index] : core_orbitals[index];
}

//******************************************************************************
int Wavefunction::maxCore_n(int ka) const
// Returns the largest n for states with kappa = ka in the core
// Note: ka is optional input; if none given, will be 0 (& not used)
// (used for energy guesses)
// Note: if you give it l instead of kappa, still works!
{
  int max_n = 0;
  for (const auto &phi : core_orbitals) {
    if (phi.k != ka && ka != 0)
      continue;
    if (phi.n > max_n)
      max_n = phi.n;
  }
  return max_n;
}
//******************************************************************************
int Wavefunction::maxCore_l() const {
  int max_l = 0;
  for (const auto &phi : core_orbitals) {
    if (phi.l() > max_l)
      max_l = phi.l();
  }
  return max_l;
}

//******************************************************************************
void Wavefunction::solveInitialCore(const std::string &str_core,
                                    int log_dele_or)
// Solves the Dirac eqn for each state in the core
// Only for local potential (direct part)
// HF/HartreeFockClass.cpp has routines for Hartree Fock
{
  if (!core_orbitals.empty()) {
    std::cerr << "WARNING 254 in Wavefunction:solveInitialCore: States "
                 "already exist! "
              << core_orbitals.size() << "\n";
    return;
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

  // Calculate c_ab = <a|b>  [only for b>a -- symmetric]
  std::vector<std::vector<double>> c_ab(Ns, std::vector<double>(Ns));
  for (std::size_t a = 0; a < Ns; a++) {
    const auto &phi_a = in_orbs[a];
    for (auto b = a + 1; b < Ns; b++) {
      const auto &phi_b = in_orbs[b];
      if (phi_a.k != phi_b.k) //|| phi_a.n == phi_b.n - can't happen!
        continue;
      c_ab[a][b] = 0.5 * (phi_a * phi_b);
    }
  }
  // note: above loop executes psia*psib half as many times as below would

  // Orthogonalise + re-norm orbitals:
  for (std::size_t a = 0; a < Ns; a++) {
    auto &phi_a = in_orbs[a];
    for (std::size_t b = 0; b < Ns; b++) {
      const auto &phi_b = in_orbs[b];
      if (phi_a.k != phi_b.k || phi_a.n == phi_b.n)
        continue;
      double cab = (a < b) ? c_ab[a][b] : c_ab[b][a];
      phi_a -= cab * phi_b;
    }
    phi_a.normalise();
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
    psi_v -= (psi_v * psi_c) * psi_c;
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
    if (l == config.l && n == config.n) {
      num_el_this = config.num;
      break;
    }
    num_el_below += config.num;
  }

  // effective Z (for energy guess) -- not perfect!
  double Zeff = 1.0 + (m_Z - num_el_below - 0.5 * num_el_this);
  if (Zeff < 1.0) {
    Zeff = 1.;
  }

  double en_a = -0.5 * std::pow(Zeff / n, 2);
  if (n > 1) {
    en_a *= 0.5;
  }
  if (Zeff < 10) {
    if (l == 0)
      en_a *= 2.5;
    if (l == 1)
      en_a *= 3.5;
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
    neff += 2. * std::pow(x, 0.5);
  if (l >= 3)
    neff += 4. * x;
  return -0.5 / std::pow(neff, 2);
}

//******************************************************************************
std::string Wavefunction::nuclearParams() const {
  std::ostringstream output;

  auto rrms = m_nuc_params.r_rms;
  auto t = m_nuc_params.t;

  switch (m_nuc_params.type) {
  case Nuclear::Type::zero:
    output << "Zero-size nucleus; ";
    break;
  case Nuclear::Type::spherical:
    output << "Spherical nucleus; "
           << " r_rms = " << rrms
           << ", r_charge = " << Nuclear::c_hdr_formula_rrms_t(rrms, 0);
    break;
  case Nuclear::Type::Fermi:
    output << "Fermi nucleus; "
           << " r_rms = " << rrms
           << ", c_hdr = " << Nuclear::c_hdr_formula_rrms_t(rrms, t)
           << ", t = " << t;
    break;
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
  if (Ncore() < 1)
    return;

  std::cout
      << "     state   k   Rinf its   eps         En (au)        En (/cm)\n";
  auto index_list = sortedEnergyList(core_orbitals, sorted);
  for (auto i : index_list) {
    const auto &phi = core_orbitals[i];
    auto r_inf = rinf(phi);
    printf("%2i) %7s %2i  %5.1f %2i  %5.0e %15.9f %15.3f", int(i),
           phi.symbol().c_str(), phi.k, r_inf, phi.its, phi.eps, phi.en,
           phi.en *PhysConst::Hartree_invcm);
    if (phi.occ_frac < 1.0)
      printf("     [%4.2f]\n", phi.occ_frac);
    else
      std::cout << "\n";
  }
  if (m_pHF) {
    auto core_energy = coreEnergyHF();
    printf("E_core = %.8g au; = %.8g /cm\n", core_energy,
           core_energy * PhysConst::Hartree_invcm);
  }
}

//******************************************************************************
void Wavefunction::printValence(
    bool sorted, const std::vector<DiracSpinor> &in_orbitals) const {
  auto tmp_orbs = (in_orbitals.empty()) ? valence_orbitals : in_orbitals;
  if (tmp_orbs.empty())
    return;

  // Find lowest valence energy:
  auto e0 = 0.0;
  for (auto &phi : tmp_orbs) {
    if (phi.en < e0)
      e0 = phi.en;
  }

  std::cout
      << "Val: state   "
      << "k   Rinf its   eps         En (au)        En (/cm)   En (/cm)\n";
  auto index_list = sortedEnergyList(tmp_orbs, sorted);
  for (auto i : index_list) {
    const auto &phi = tmp_orbs[i];
    auto r_inf = rinf(phi);
    printf("%2i) %7s %2i  %5.1f %2i  %5.0e %15.9f %15.3f", int(i),
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

//******************************************************************************
std::vector<double> Wavefunction::coreDensity() const {
  std::vector<double> rho(rgrid.ngp, 0.0);
  for (const auto &phi : core_orbitals) {
    auto f = double(phi.twoj() + 1) * phi.occ_frac;
    for (auto i = 0ul; i < rgrid.ngp; i++) {
      rho[i] += f * (phi.f[i] * phi.f[i] + phi.g[i] * phi.g[i]);
    }
  }
  return rho;
}
