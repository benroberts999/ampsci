#include "Wavefunction/Wavefunction.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp" //just for enum..
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/FeynmanSigma.hpp"
#include "MBPT/GoldstoneSigma.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/Parametric_potentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

//******************************************************************************
Wavefunction::Wavefunction(const GridParameters &gridparams,
                           const Nuclear::Parameters &nuc_params,
                           double var_alpha)
    : rgrid(std::make_shared<const Grid>(gridparams)),
      alpha(PhysConst::alpha * var_alpha),
      m_nuclear(nuc_params),
      vnuc(Nuclear::formPotential(nuc_params, rgrid->r())) {
  if (alpha * m_nuclear.z > 1.0) {
    std::cerr << "Alpha too large: Z*alpha=" << m_nuclear.z * alpha << "\n";
    std::abort();
  }
}

//******************************************************************************
Wavefunction::Wavefunction(const Wavefunction &wf)
    : Wavefunction(wf.rgrid->params(), wf.get_nuclearParameters(),
                   wf.alpha / PhysConst::alpha) {
  // NOTE: orbitals in new_wf point to OLD grid (*(wf.rgrid), not this->rgrid)
  // new WF ONLY has orbitals, does not have HF/Sigma etc!
  this->core = wf.core;
  this->valence = wf.valence;
  this->basis = wf.basis;
  this->spectrum = wf.spectrum;
  this->vnuc = wf.vnuc;
  this->vdir = wf.vdir;
  if (wf.qed)
    this->qed = std::make_unique<QED::RadPot>(*wf.qed);
  this->m_core_configs = wf.m_core_configs;
  this->num_core_electrons = wf.num_core_electrons;
  this->m_core_string = wf.m_core_string;
}

//******************************************************************************
void Wavefunction::solveDirac(DiracSpinor &psi, double e_a,
                              const std::vector<double> &vex,
                              int log_dele_or) const
// Uses Adams::boundState to solve Dirac Eqn for local potential (Vnuc + Vdir)
// If no e_a is given, will use the existing one!
// (Usually, a better guess should be given, using P.T.)
// Note: optionally takes in exchange potential! (see overloaded above)
// Note: Uses the "dodgy" re-scaled exchange potenital:
// Vex\psi_a = sum_b vex_a psi_b -> [sum_b vex_a (psi_b/psi_a)] psi_a
// so, here, vex = [sum_b vex_a (psi_b/psi_a)]
// This is not ideal..
{
  const auto v_a = qip::add(get_Vlocal(psi.l()), vex);
  if (e_a != 0) {
    psi.set_en() = e_a;
  } else if (psi.en() == 0) {
    psi.set_en() = enGuessVal(psi.n, psi.k);
  }
  DiracODE::boundState(psi, psi.en(), v_a, get_Hmag(psi.l()), alpha,
                       log_dele_or);
}

//------------------------------------------------------------------------------
void Wavefunction::solveDirac(DiracSpinor &psi, double e_a,
                              int log_dele_or) const
// Overloaded version; see above
// This one doesn't have exchange potential
{
  // std::vector<double> empty_vec;
  return solveDirac(psi, e_a, {}, log_dele_or);
}

//******************************************************************************
void Wavefunction::determineCore(const std::string &str_core_in)
// Takes in a string list for the core configuration, outputs an int list
// Takes in previous closed shell (noble), + 'rest' (or just the rest)
// E.g:
//   Core of Cs: Xe
//   Core of Gold: Xe 4f14 5d10
// 'rest' is in form nLm : n=n, L=l, m=number of electrons in that nl shell.
{
  m_core_configs = AtomData::core_parser(str_core_in);

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

  if (num_core_electrons > m_nuclear.z) {
    std::cout << "Problem with core: " << str_core_in << "\n";
    std::cout << "= " << m_core_string << "\n";
    std::cout << "= " << AtomData::niceCoreOutput(m_core_string) << "\n";
    std::cout << "Too many electrons: N_core=" << num_core_electrons
              << ", Z=" << m_nuclear.z << "\n";
    std::abort();
  }

  return;
}

//******************************************************************************
// void Wavefunction::solve_core(const std::string &method,
//                                    const double x_Breit,
//                                    const std::string &in_core, double eps_HF,
//                                    bool print) {
//   if (m_pHF == nullptr) {
//     solveLocalCore(in_core, 16);
//     m_pHF = std::make_unique<HF::HartreeFock>(this, HF::parseMethod(method),
//                                               x_Breit, eps_HF);
//   }
//   m_pHF->verbose = print;
//   vdir = m_pHF->solve_core();
// }

void Wavefunction::solve_core(const std::string &method, const double x_Breit,
                              const std::string &in_core, double eps_HF,
                              bool print) {
  // core config
  // initial core if HF

  if (!m_pHF) {

    if (in_core != "") {
      double h_g = 0, d_t = 0;
      // XXX Update to use other coefs?
      if (method != "Local") {
        Parametric::defaultGreenCore(m_nuclear.z, h_g, d_t);
      } else {
        Parametric::defaultGreen(m_nuclear.z, h_g, d_t);
      }
      vdir = Parametric::GreenPotential(m_nuclear.z, rgrid->r(), h_g, d_t);
    }
    const auto eps_targ = method == "Local" ? 16 : 5;
    solveLocalCore(in_core, eps_targ);
    m_pHF = std::make_unique<HF::HartreeFock>(this, HF::parseMethod(method),
                                              x_Breit, eps_HF);
    m_pHF->verbose = print;
  }

  if (method != "Local")
    vdir = m_pHF->solve_core();
}

//******************************************************************************
auto Wavefunction::coreEnergyHF() const {
  if (!m_pHF) {
    return 0.0;
  }
  return m_pHF->calculateCoreEnergy();
}

//******************************************************************************
void Wavefunction::solve_valence(const std::string &in_valence_str,
                                 const bool print) {
  if (!m_pHF) {
    std::cerr << "WARNING 62: Cant call solve_valence before "
                 "solve_core\n";
    return;
  }

  auto eps = m_pHF->method() == HF::Method::Local ? 16 : 5;

  if (m_pHF->method() == HF::Method::KohnSham) {
    localValence(in_valence_str, true);
  } else {
    const auto val_lst = AtomData::listOfStates_nk(in_valence_str);
    for (const auto &[n, k, en] : val_lst) {
      (void)en;
      if (!isInCore(n, k) && !isInValence(n, k)) {
        solveNewValence(n, k, 0, eps);
      }
    }
    if (m_pHF->method() != HF::Method::Local)
      m_pHF->solveValence(&valence, print);
  }
}

//******************************************************************************
void Wavefunction::localValence(const std::string &in_valence_str,
                                bool list_each) {

  // Use for Kohn-Sham.
  auto val_lst = list_each ? AtomData::listOfStates_singlen(in_valence_str) :
                             AtomData::listOfStates_nk(in_valence_str);
  for (const auto &[n, k, en] : val_lst) {
    (void)en;
    solveNewValence(n, k, 0, 17);
  }
}
//******************************************************************************
void Wavefunction::radiativePotential(QED::RadPot::Scale scale, double rcut,
                                      double scale_rN,
                                      const std::vector<double> &x_spd,
                                      bool do_readwrite, bool print) {

  if (x_spd.empty())
    return;

  const auto r_N_au =
      std::sqrt(5.0 / 3.0) * scale_rN * m_nuclear.r_rms / PhysConst::aB_fm;

  qed = std::make_unique<QED::RadPot>(QED::RadPot(
      rgrid->r(), Znuc(), r_N_au, rcut, scale, x_spd, print, do_readwrite));

  // If HF already exists, update it to include new qed!
  if (m_pHF)
    m_pHF->update_Vrad(qed.get());
}

//******************************************************************************
bool Wavefunction::isInCore(int n, int k) const {
  const auto find_nk = [n, k](const auto &Fa) {
    return Fa.n == n && Fa.k == k;
  };
  const auto Fnk = std::find_if(cbegin(core), cend(core), find_nk);
  return Fnk != cend(core);
}
//------------------------------------------------------------------------------
bool Wavefunction::isInValence(int n, int k) const {
  const auto find_nk = [n, k](const auto Fa) { return Fa.n == n && Fa.k == k; };
  const auto Fnk = std::find_if(cbegin(valence), cend(valence), find_nk);
  return Fnk != cend(valence);
}

//******************************************************************************
const DiracSpinor *Wavefunction::getState(int n, int k,
                                          bool *is_valence) const {
  const auto find_nk = [n, k](const auto Fa) { return Fa.n == n && Fa.k == k; };
  // Try to find in core:
  auto Fnk = std::find_if(cbegin(core), cend(core), find_nk);
  if (Fnk != cend(core)) {
    if (is_valence != nullptr)
      *is_valence = false;
    return &*Fnk;
  }
  // If not in core, try to find in valence:
  Fnk = std::find_if(cbegin(valence), cend(valence), find_nk);
  if (Fnk != cend(valence)) {
    if (is_valence != nullptr)
      *is_valence = true;
    return &*Fnk;
  }
  // otherwise, return nope
  return nullptr;
}

const DiracSpinor *Wavefunction::getState(std::string_view state,
                                          bool *is_valence) const {
  // std::pair<int, int> parse_symbol(std::string_view symbol);
  const auto [n, k] = AtomData::parse_symbol(state);
  return getState(n, k, is_valence);
}

//******************************************************************************
int Wavefunction::maxCore_n(int ka) const
// Returns the largest n for states with kappa = ka in the core
// Note: ka is optional input; if none given, will be 0 (& not used)
// (used for energy guesses)
// Note: if you give it l instead of kappa, still works!
{
  int max_n = 0;
  for (const auto &phi : core) {
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
  for (const auto &phi : core) {
    if (phi.l() > max_l)
      max_l = phi.l();
  }
  return max_l;
}

//******************************************************************************
double Wavefunction::en_coreval_gap() const {
  // Find core/valence energy: allows distingush core/valence states
  const auto ec_max = core.empty() ? 0.0 :
                                     std::max_element(cbegin(core), cend(core),
                                                      DiracSpinor::comp_en)
                                         ->en();
  const auto ev_min =
      valence.empty() ?
          0.0 :
          std::min_element(cbegin(valence), cend(valence), DiracSpinor::comp_en)
              ->en();
  return 0.5 * (ev_min + ec_max);
}

//******************************************************************************
void Wavefunction::solveLocalCore(const std::string &str_core, int log_dele_or)
// Solves the Dirac eqn for each state in the core
// Only for local potential (direct part)
// HF/HartreeFock.cpp has routines for Hartree Fock
{
  if (!core.empty()) {
    core.clear();           //?
    m_core_configs.clear(); //?
  }

  determineCore(str_core); // sets m_core_configs :( ?

  for (const auto &[n, l, num] : m_core_configs) {
    if (num == 0)
      continue;
    double en_a = enGuessCore(n, l);
    int k1 = l; // j = l-1/2
    if (k1 != 0) {
      auto &new_Fc = core.emplace_back(n, k1, rgrid);
      solveDirac(new_Fc, en_a, log_dele_or);
      new_Fc.set_occ_frac() = double(num) / (4 * l + 2);
      en_a = 0.95 * new_Fc.en();
      if (en_a > 0)
        en_a = enGuessCore(n, l);
    }
    int k2 = -(l + 1); // j=l+1/2
    auto &new_Fc = core.emplace_back(n, k2, rgrid);
    solveDirac(new_Fc, en_a, log_dele_or);
    new_Fc.set_occ_frac() = double(num) / (4 * l + 2);
  }
}

//******************************************************************************
void Wavefunction::solveNewValence(int n, int k, double en_a, int log_dele_or)
// Update to take a list ok nken's ?
{
  valence.emplace_back(n, k, rgrid);

  // Solve local dirac Eq:
  auto &psi = valence.back();
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
void Wavefunction::orthogonaliseWrt(DiracSpinor &psi_v,
                                    const std::vector<DiracSpinor> &in_orbs) {
  for (const auto &psi_c : in_orbs) {
    if (psi_v.k != psi_c.k)
      continue;
    psi_v -= (psi_v * psi_c) * psi_c;
  }
}
//******************************************************************************
void Wavefunction::orthonormaliseWrt(DiracSpinor &psi_v,
                                     const std::vector<DiracSpinor> &in_orbs)
// Static.
// Force given orbital to be orthogonal to all core orbitals
// [After the core is 'frozen', don't touch core orbitals!]
// |v> --> |v> - sum_c |c><c|v>
// note: here, c denotes core orbitals
{
  orthogonaliseWrt(psi_v, in_orbs);
  psi_v.normalise();
}

//******************************************************************************
std::tuple<double, double> Wavefunction::lminmax_core_range(int l,
                                                            double eps) const {
  std::vector<double> rho_l(rgrid->num_points());
  bool found = false;
  for (const auto &Fc : core) {
    if (Fc.l() == l || l < 0) {
      found = true;
      qip::add(&rho_l, Fc.rho());
    }
  }
  if (!found) {
    return {rgrid->r().front(), rgrid->r().back()};
  }
  // find maximum rho:
  const auto max = std::max_element(rho_l.begin(), rho_l.end());
  // find first position that rho=|psi^2| reaches cut-off
  const auto cut = max != rho_l.end() ? eps * (*max) : 0.0;
  const auto lam = [cut](const auto &v) { return v >= cut; };
  const auto first = std::find_if(rho_l.begin(), rho_l.end(), lam);
  const auto last = std::find_if(rho_l.rbegin(), rho_l.rend(), lam);
  const auto index_first = std::size_t(first - rho_l.begin());
  const auto index_last = std::size_t(rho_l.rend() - last);
  return {rgrid->r()[index_first], rgrid->r()[index_last]};
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
  double Zeff = 1.0 + (m_nuclear.z - num_el_below - 0.5 * num_el_this);
  if (Zeff < 1.0) {
    Zeff = 1.0;
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
  const int maxn = maxCore_n();
  const int l = AtomData::l_k(ka);
  const int dn = n - maxn;
  double neff = 1.0 + dn;
  double x = 1;
  double Z_eff = m_nuclear.z - num_core_electrons;
  if (Z_eff <= 0)
    Z_eff = 0.5;
  if (maxn < 4)
    x = 0.25;
  if (l == 1)
    neff += 0.5 * x;
  if (l == 2)
    neff += 2.0 * std::pow(x, 0.5);
  if (l >= 3)
    neff += 4.0 * x;
  return -0.5 * Z_eff * Z_eff / std::pow(neff, 2);
}

//******************************************************************************
std::string Wavefunction::nuclearParams() const {
  std::ostringstream output;

  auto rrms = m_nuclear.r_rms;
  auto t = m_nuclear.t;

  switch (m_nuclear.type) {
  case Nuclear::Type::point:
    output << "Point-like nucleus; ";
    break;
  case Nuclear::Type::spherical:
    output << "Spherical nucleus; "
           << " r_rms = " << rrms
           << ", r_charge = " << Nuclear::c_hdr_formula_rrms_t(rrms, 0);
    break;
  case Nuclear::Type::Gaussian:
    output << "Gaussian nucleus; "
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
                               bool do_sort)
// Static
// Outouts a list of integers corresponding to the states
// sorted by energy (lowest energy first)
{
  using DoubleInt = std::pair<double, std::size_t>;
  std::vector<DoubleInt> t_en;
  for (std::size_t i = 0; i < tmp_orbs.size(); i++) {
    t_en.emplace_back(tmp_orbs[i].en(), i);
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

  if (Ncore() == 0) {
    std::cout << " - H-like";
  }
  std::cout << ")\n";
  if (Ncore() < 1)
    return;

  std::cout
      << "     state   k   Rinf its   eps         En (au)        En (/cm)\n";
  auto index_list = sortedEnergyList(core, sorted);
  for (auto i : index_list) {
    const auto &phi = core[i];
    auto r_inf = rgrid->r()[phi.max_pt() - 1]; // rinf(phi);
    printf("%2i) %7s %2i  %5.1f %2i  %5.0e %15.9f %15.3f", int(i),
           phi.symbol().c_str(), phi.k, r_inf, phi.its(), phi.eps(), phi.en(),
           phi.en() * PhysConst::Hartree_invcm);
    if (phi.occ_frac() < 1.0)
      printf("     [%4.2f]\n", phi.occ_frac());
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
  auto tmp_orbs = (in_orbitals.empty()) ? valence : in_orbitals;
  if (tmp_orbs.empty())
    return;

  // Find lowest valence energy:
  auto e0 = 0.0;
  for (auto &phi : tmp_orbs) {
    if (phi.en() < e0)
      e0 = phi.en();
  }

  std::cout
      << "Val: state   "
      << "k   Rinf its   eps         En (au)        En (/cm)   En (/cm)\n";
  auto index_list = sortedEnergyList(tmp_orbs, sorted);
  for (auto i : index_list) {
    const auto &phi = tmp_orbs[i];
    auto r_inf = rgrid->r()[phi.max_pt() - 1]; // rinf(phi);
    printf("%2i) %7s %2i  %5.1f %2i  %5.0e %15.9f %15.3f", int(i),
           phi.symbol().c_str(), phi.k, r_inf, phi.its(), phi.eps(), phi.en(),
           phi.en() * PhysConst::Hartree_invcm);
    printf(" %10.2f\n", (phi.en() - e0) * PhysConst::Hartree_invcm);
  }
}

//******************************************************************************
void Wavefunction::printBasis(const std::vector<DiracSpinor> &the_basis,
                              bool sorted) const {
  std::cout
      << "     State   k  R0     Rinf          En(Basis)         En(HF)\n";
  const auto index_list = sortedEnergyList(the_basis, sorted);
  for (const auto i : index_list) {
    const auto &phi = the_basis[i];
    const auto r_0 = phi.r0();
    const auto r_inf = phi.rinf();

    const auto *hf_phi = getState(phi.n, phi.k);
    if (hf_phi != nullptr) {
      // found HF state
      const auto eps =
          2.0 * (phi.en() - hf_phi->en()) / (phi.en() + hf_phi->en());
      printf("%2i) %7s %2i  %5.0e %5.1f %18.7f  %13.7f  %6.0e\n", int(i),
             phi.symbol().c_str(), phi.k, r_0, r_inf, phi.en(), hf_phi->en(),
             eps);
    } else {
      printf("%2i) %7s %2i  %5.0e %5.1f %18.7f\n", int(i), phi.symbol().c_str(),
             phi.k, r_0, r_inf, phi.en());
    }
  }
}

//******************************************************************************
std::vector<double> Wavefunction::coreDensity() const {
  std::vector<double> rho(rgrid->num_points(), 0.0);
  for (const auto &phi : core) {
    auto f = double(phi.twoj() + 1) * phi.occ_frac();
    for (auto i = 0ul; i < rgrid->num_points(); i++) {
      rho[i] += f * (phi.f(i) * phi.f(i) + phi.g(i) * phi.g(i));
    }
  }
  return rho;
}

//******************************************************************************
void Wavefunction::formBasis(const SplineBasis::Parameters &params) {
  if (params.n > 0) {
    IO::ChronoTimer t("Basis");
    basis = SplineBasis::form_basis(params, *this, false);
  }
}
//------------------------------------------------------------------------------
void Wavefunction::formSpectrum(const SplineBasis::Parameters &params) {
  if (params.n > 0) {
    IO::ChronoTimer t("Spectrum");
    spectrum = SplineBasis::form_basis(params, *this, true);
  }
}

//******************************************************************************
void Wavefunction::formSigma(
    const int nmin_core, const bool form_matrix, const double r0,
    const double rmax, const int stride, const bool each_valence,
    const bool include_G, const std::vector<double> &lambdas,
    const std::vector<double> &fk, const std::string &in_fname,
    const std::string &out_fname, const bool FeynmanQ, const bool ScreeningQ,
    const bool holeParticleQ, const int lmax, const bool GreenBasis,
    const bool PolBasis, const double omre, double w0, double wratio,
    const std::optional<IO::InputBlock> &ek) {
  if (valence.empty())
    return;

  /*
      // XXX Re-factor
      a) Make sub-block for Feynman, Goldstone Options
      b) make sub-block for fit_to: e.g., fitTo = [6s+=31406;];
  */
  const std::string ext = FeynmanQ ? ".sigf" : ".sig2";
  const auto ifname = in_fname == "" ? identity() + ext : in_fname + ext;
  const auto ofname = out_fname == "" ? identity() + ext : out_fname + ext;

  const auto method =
      FeynmanQ ? MBPT::Method::Feynman : MBPT::Method::Goldstone;

  const auto sigp = MBPT::Sigma_params{
      method, nmin_core, include_G, lmax,       GreenBasis,    PolBasis,
      omre,   w0,        wratio,    ScreeningQ, holeParticleQ, fk};

  const auto subgridp = MBPT::rgrid_params{r0, rmax, std::size_t(stride)};

  // Correlaion potential matrix:
  switch (method) {
  case MBPT::Method::Goldstone:
    m_Sigma = std::make_unique<MBPT::GoldstoneSigma>(m_pHF.get(), basis, sigp,
                                                     subgridp, ifname);
    break;
  case MBPT::Method::Feynman:
    m_Sigma = std::make_unique<MBPT::FeynmanSigma>(m_pHF.get(), basis, sigp,
                                                   subgridp, ifname);
    break;
  }

  // This is for each valence state.... otherwise, just do for lowest??
  if (form_matrix && !valence.empty()) {
    if (each_valence) {
      // calculate sigma for each valence state:
      for (const auto &Fv : valence) {
        m_Sigma->formSigma(Fv.k, Fv.en(), Fv.n);
      }
    } else if (!ek && m_Sigma->empty()) {
      // calculate sigma for lowest n valence state of each kappa:
      const auto max_ki = DiracSpinor::max_kindex(valence);
      for (int ki = 0; ki <= max_ki; ++ki) {
        auto Fv = std::find_if(cbegin(valence), cend(valence),
                               [ki](auto f) { return f.k_index() == ki; });
        if (Fv != cend(valence))
          m_Sigma->formSigma(Fv->k, Fv->en(), Fv->n);
      }
    } else if (ek && m_Sigma->empty()) {
      // solve at specific energies:
      for (auto &[state, en] : ek->options()) {
        auto [n, k] = AtomData::parse_symbol(state);
        m_Sigma->formSigma(k, std::stod(en), n);
      }
    }
  }

  if (!lambdas.empty()) {
    std::cout << "Note: be careful order of input lambdas matches Sigma\n";
    // If Sigma diverged from valence, may be diff order to valence...
    // Better to make it explicit!
    m_Sigma->scale_Sigma(lambdas);
    m_Sigma->print_scaling();
  }

  if (out_fname != "false" && form_matrix)
    m_Sigma->read_write(ofname, IO::FRW::RoW::write);
}

//******************************************************************************
void Wavefunction::hartreeFockBrueckner(const bool print) {
  if (!m_pHF) {
    std::cerr << "WARNING 62: Cant call solve_valence before "
                 "solve_core\n";
    return;
  }
  if (m_Sigma)
    m_pHF->solveBrueckner(&valence, *(m_Sigma.get()), print);
}

//******************************************************************************
void Wavefunction::fitSigma_hfBrueckner(
    const std::string &, const std::vector<double> &fit_energies) {
  std::cout << "\nFitting Sigma for lowest valence states:\n";

  const auto max_its = 10;
  const auto eps_targ = 1.0e-7;

  // XXX Assume the 'fit_to' are in same order as valence!!
  // Must be called after HF, before Bruckenr....

  //
  for (auto i = 0ul; i < fit_energies.size(); ++i) {
    if (i >= valence.size())
      break;
    const auto &Fv = valence[i];
    const auto e_exp = fit_energies[i];
    if (e_exp >= 0.0)
      continue;

    // std::cout << Fv.symbol() << " " << Fv.en() * PhysConst::Hartree_invcm
    //           << " -> " << e_exp * PhysConst::Hartree_invcm << ": ";

    printf("%4s %7.0f [%7.0f] : ", Fv.shortSymbol().c_str(),
           Fv.en() * PhysConst::Hartree_invcm,
           e_exp * PhysConst::Hartree_invcm);
    const double en_0 = Fv.en(); // HF value
    auto lambda = 1.0;
    double eps = 1.0;
    int its = 0;
    for (; its < max_its; its++) {
      auto Fv_l = Fv;
      m_Sigma->scale_Sigma(Fv_l.n, Fv_l.k, lambda);
      // nb: hf_Brueckner must start from HF... so, call on copy of Fv....
      m_pHF->hf_Brueckner(Fv_l, *m_Sigma);
      double en_l = Fv_l.en();
      eps = std::abs((e_exp - en_l) / e_exp);
      if (eps < eps_targ)
        break;
      const auto r = (e_exp - en_l) / (en_l - en_0);
      lambda = std::clamp(lambda * (1.0 + r), 0.5, 1.5);
    }
    printf("%.0e (%2i); lambda = %.4f\n", eps, its, lambda);
  }

  hartreeFockBrueckner(true);
}

//******************************************************************************
void Wavefunction::SOEnergyShift() {
  std::cout << "\nMBPT(2): Second-order valence energy shifts\n";
  std::cout << "and matrix elements <v|Sigma(2)|v>:\n";

  if (!m_Sigma) {
    std::cout << "No Sigma?\n";
    return;
  }

  double e0 = 0;
  std::cout << "state |  E(HF)      E(2)       <v|S2|v> |  E(HF+2)     E(HF+2) "
               " (cm^-1)\n";
  for (const auto &v : valence) {
    const auto delta = m_Sigma->SOEnergyShift(v, v);
    const auto delta2 = v * (*m_Sigma)(v);
    const auto cm = PhysConst::Hartree_invcm;
    if (e0 == 0)
      e0 = (v.en() + delta);
    printf("%6s| %9.6f  %+9.6f  %+9.6f | %9.6f = %8.1f  %7.1f\n",
           v.symbol().c_str(), v.en(), delta, delta2, (v.en() + delta),
           (v.en() + delta) * cm, (v.en() + delta - e0) * cm);
    if (std::abs(delta / v.en()) > 0.2)
      std::cout << "      *** Warning: delta too large?\n";
  }
}

//******************************************************************************
std::vector<double> Wavefunction::get_Vlocal(int l) const {
  return qed ? qip::add(vnuc, vdir, qed->Vel(l)) : qip::add(vnuc, vdir);
}
//******************************************************************************
std::vector<double> Wavefunction::get_Hmag(int l) const {
  return qed ? qed->Hmag(l) : std::vector<double>{};
}

//******************************************************************************
double Wavefunction::Hab(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
  if (Fa.k != Fb.k)
    return 0.0;
  const auto kappa = Fa.k;
  const auto max = std::min(Fa.max_pt(), Fb.max_pt());
  const auto min = std::max(Fa.min_pt(), Fb.min_pt());
  const auto &drdu = Fa.rgrid->drdu();

  const auto the_same = &Fa == &Fb;

  auto dga = NumCalc::derivative(Fa.g(), drdu, Fb.rgrid->du(), 1);
  auto dgb =
      the_same ? dga : NumCalc::derivative(Fb.g(), drdu, Fb.rgrid->du(), 1);

  for (std::size_t i = min; i < max; i++) {
    const auto r = Fa.rgrid->r(i);
    dga[i] -= (kappa * Fa.g(i) / r);
    dgb[i] -= (kappa * Fb.g(i) / r);
  }

  const auto D1m2 = NumCalc::integrate(1.0, min, max, Fa.f(), dgb, drdu) +
                    NumCalc::integrate(1.0, min, max, Fb.f(), dga, drdu);

  const auto Sab = NumCalc::integrate(1.0, min, max, Fa.g(), Fb.g(), drdu);

  const auto &v = get_Vlocal(Fa.l());
  const auto Vab = NumCalc::integrate(1.0, min, max, Fa.f(), Fb.f(), v, drdu) +
                   NumCalc::integrate(1.0, min, max, Fa.g(), Fb.g(), v, drdu);

  const auto &Hmag = get_Hmag(Fa.l());

  const auto V_mag =
      Hmag.empty() ?
          0.0 :
          NumCalc::integrate(1.0, min, max, Fa.f(), Fb.g(), Hmag, drdu) +
              NumCalc::integrate(1.0, min, max, Fa.g(), Fb.f(), Hmag, drdu);
  const auto c = 1.0 / alpha;

  return (Vab - c * (D1m2 + 2.0 * c * Sab + V_mag)) * Fa.rgrid->du();
}

double Wavefunction::Hab(const DiracSpinor &Fa, const DiracSpinor &dFa,
                         const DiracSpinor &Fb, const DiracSpinor &dFb) const {
  if (Fa.k != Fb.k)
    return 0.0;
  const auto kappa = Fa.k;
  const auto max = std::min(Fa.max_pt(), Fb.max_pt());
  const auto min = std::max(Fa.min_pt(), Fb.min_pt());
  const auto &drdu = Fa.rgrid->drdu();

  auto dga = dFa.g();
  auto dgb = dFb.g();

  for (std::size_t i = min; i < max; i++) {
    const auto r = Fa.rgrid->r(i);
    dga[i] -= (kappa * Fa.g(i) / r);
    dgb[i] -= (kappa * Fb.g(i) / r);
  }

  const auto D1m2 = NumCalc::integrate(1.0, min, max, Fa.f(), dgb, drdu) +
                    NumCalc::integrate(1.0, min, max, Fb.f(), dga, drdu);

  const auto Sab = NumCalc::integrate(1.0, min, max, Fa.g(), Fb.g(), drdu);

  const auto &v = get_Vlocal(Fa.l());
  const auto Vab = NumCalc::integrate(1.0, min, max, Fa.f(), Fb.f(), v, drdu) +
                   NumCalc::integrate(1.0, min, max, Fa.g(), Fb.g(), v, drdu);

  const auto &Hmag = get_Hmag(Fa.l());

  const auto V_mag =
      Hmag.empty() ?
          0.0 :
          NumCalc::integrate(1.0, min, max, Fa.f(), Fb.g(), Hmag, drdu) +
              NumCalc::integrate(1.0, min, max, Fa.g(), Fb.f(), Hmag, drdu);
  const auto c = 1.0 / alpha;

  return (Vab - c * (D1m2 + 2.0 * c * Sab + V_mag)) * Fa.rgrid->du();
}
