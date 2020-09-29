#include "Wavefunction/Wavefunction.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp" //just for enum..
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/FeynmanSigma.hpp"
#include "MBPT/GoldstoneSigma2.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/Parametric_potentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/RadiativePotential.hpp"
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
      vnuc(Nuclear::formPotential(nuc_params, rgrid->r)) {
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
  this->vrad = wf.vrad;
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
    psi.en = e_a;
  } else if (psi.en == 0) {
    psi.en = enGuessVal(psi.n, psi.k);
  }
  DiracODE::boundState(psi, psi.en, v_a, get_Hmag(psi.l()), alpha, log_dele_or);
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
void Wavefunction::hartreeFockCore(const std::string &method,
                                   const double x_Breit,
                                   const std::string &in_core, double eps_HF) {
  solveInitialCore(in_core, 5);
  m_pHF = std::make_unique<HF::HartreeFock>(this, HF::parseMethod(method),
                                            x_Breit, eps_HF);
  vdir = m_pHF->solveCore();
}

//******************************************************************************
auto Wavefunction::coreEnergyHF() const {
  if (!m_pHF) {
    return 0.0;
  }
  return m_pHF->calculateCoreEnergy();
}

//******************************************************************************
void Wavefunction::hartreeFockValence(const std::string &in_valence_str,
                                      const bool print) {
  if (!m_pHF) {
    std::cerr << "WARNING 62: Cant call hartreeFockValence before "
                 "hartreeFockCore\n";
    return;
  }
  auto val_lst = AtomData::listOfStates_nk(in_valence_str);
  for (const auto &[n, k, en] : val_lst) {
    (void)en;
    if (!isInCore(n, k) && !isInValence(n, k)) {
      solveNewValence(n, k, 0, 5);
    }
  }
  m_pHF->solveValence(&valence, print);
}

//******************************************************************************
void Wavefunction::localValence(const std::string &in_valence_str,
                                bool list_each) {

  // Use for Kohn-Sham.
  auto val_lst = list_each ? AtomData::listOfStates_singlen(in_valence_str)
                           : AtomData::listOfStates_nk(in_valence_str);
  for (const auto &[n, k, en] : val_lst) {
    (void)en;
    solveNewValence(n, k, 0, 17);
  }
}
//******************************************************************************
void Wavefunction::radiativePotential(double x_simple, double x_Ueh,
                                      double x_SEe_h, double x_SEe_l,
                                      double x_SEm, double rcut,
                                      double scale_rN,
                                      const std::vector<double> &x_spd) {
  if (x_spd.empty())
    return;

  const auto r_rms_fm = scale_rN * m_nuclear.r_rms;

  auto is_one = [](double x) { return std::abs(x - 1.0) < 1.0e-6; };
  // auto is_zero = [](double x) { return std::abs(x) < 1.0e-6; };

  std::string label = "_";
  if (is_one(x_Ueh))
    label += 'u';
  if (is_one(x_SEe_h))
    label += 'h';
  if (is_one(x_SEe_l))
    label += 'l';
  if (is_one(x_SEm))
    label += 'm';
  if (is_one(x_simple))
    label += "_simple";

  const auto fname = atomicSymbol() + label + ".qed";
  std::vector<double> Hel_tmp, Hmag_tmp;

  // Only read/write is "normal/full" rad pot is used
  const bool do_rw = true;

  const auto read_ok =
      do_rw ? RadiativePotential::read_write_qed(rgrid->r, Hel_tmp, Hmag_tmp,
                                                 fname, IO::FRW::read)
            : false;

  if (!read_ok) {
    if (x_Ueh > 0 || x_SEe_h > 0 || x_SEe_l > 0 || x_SEm > 0 || x_simple > 0) {
      std::cout
          << "\nIncluding Ginges/Flambaum QED radiative potential (up to r="
          << rcut << "):\n";
    } else if (std::abs(x_simple) > 0) {
      std::cout << "\nIncluding simple exponential radiative potential\n";
    } else {
      return;
    }
    Hel_tmp = RadiativePotential::form_Hel(rgrid->r, x_simple, x_Ueh, x_SEe_h,
                                           x_SEe_l, r_rms_fm, m_nuclear.z,
                                           alpha, rcut);
    Hmag_tmp = RadiativePotential::form_Hmag(rgrid->r, x_SEm, r_rms_fm,
                                             m_nuclear.z, alpha, rcut);
    if (do_rw)
      RadiativePotential::read_write_qed(rgrid->r, Hel_tmp, Hmag_tmp, fname,
                                         IO::FRW::write);
  }

  const bool print_xl = (x_spd.size() > 1 || x_spd.front() != 1.0);
  int l = 0;
  for (const auto &x_l : x_spd) {
    if (print_xl)
      std::cout << "l=" << l << ", x_l=" << x_l << "\n";
    auto V_l = Hel_tmp; // make a copy..dumb
    auto H_l = Hmag_tmp;
    qip::scale(&V_l, x_l);
    qip::scale(&H_l, x_l);
    vrad.set_Hel(V_l, l);
    vrad.set_Hmag(H_l, l);
    l++;
  }
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
void Wavefunction::solveInitialCore(const std::string &str_core,
                                    int log_dele_or)
// Solves the Dirac eqn for each state in the core
// Only for local potential (direct part)
// HF/HartreeFock.cpp has routines for Hartree Fock
{
  if (!core.empty()) {
    core.clear();           //?
    m_core_configs.clear(); //?
  }

  determineCore(str_core); // sets m_core_configs :( ?

  if (!m_core_configs.empty()) {
    double h_g = 0, d_t = 0;
    Parametric::defaultGreenCore(m_nuclear.z, h_g, d_t);
    vdir = Parametric::GreenPotential(m_nuclear.z, rgrid->r, h_g, d_t);
  }

  for (const auto &[n, l, num] : m_core_configs) {
    if (num == 0)
      continue;
    double en_a = enGuessCore(n, l);
    int k1 = l; // j = l-1/2
    if (k1 != 0) {
      auto &new_Fc = core.emplace_back(n, k1, rgrid);
      solveDirac(new_Fc, en_a, log_dele_or);
      new_Fc.occ_frac = double(num) / (4 * l + 2);
      en_a = 0.95 * new_Fc.en;
      if (en_a > 0)
        en_a = enGuessCore(n, l);
    }
    int k2 = -(l + 1); // j=l+1/2
    auto &new_Fc = core.emplace_back(n, k2, rgrid);
    solveDirac(new_Fc, en_a, log_dele_or);
    new_Fc.occ_frac = double(num) / (4 * l + 2);
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
  std::vector<double> rho_l(rgrid->num_points);
  bool found = false;
  for (const auto &Fc : core) {
    if (Fc.l() == l || l < 0) {
      found = true;
      qip::add(&rho_l, Fc.rho());
    }
  }
  if (!found) {
    return {rgrid->r.front(), rgrid->r.back()};
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
  return {rgrid->r[index_first], rgrid->r[index_last]};
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
    auto r_inf = rgrid->r[phi.pinf - 1]; // rinf(phi);
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
  auto tmp_orbs = (in_orbitals.empty()) ? valence : in_orbitals;
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
    auto r_inf = rgrid->r[phi.pinf - 1]; // rinf(phi);
    printf("%2i) %7s %2i  %5.1f %2i  %5.0e %15.9f %15.3f", int(i),
           phi.symbol().c_str(), phi.k, r_inf, phi.its, phi.eps, phi.en,
           phi.en *PhysConst::Hartree_invcm);
    printf(" %10.2f\n", (phi.en - e0) * PhysConst::Hartree_invcm);
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
    const auto r_0 = rgrid->r[phi.p0];
    const auto r_inf = rgrid->r[phi.pinf - 1];

    const auto *hf_phi = getState(phi.n, phi.k);
    if (hf_phi != nullptr) {
      // found HF state
      const auto eps = 2.0 * (phi.en - hf_phi->en) / (phi.en + hf_phi->en);
      printf("%2i) %7s %2i  %5.0e %5.1f %18.7f  %13.7f  %6.0e\n", int(i),
             phi.symbol().c_str(), phi.k, r_0, r_inf, phi.en, hf_phi->en, eps);
    } else {
      printf("%2i) %7s %2i  %5.0e %5.1f %18.7f\n", int(i), phi.symbol().c_str(),
             phi.k, r_0, r_inf, phi.en);
    }
  }
}

//******************************************************************************
std::vector<double> Wavefunction::coreDensity() const {
  std::vector<double> rho(rgrid->num_points, 0.0);
  for (const auto &phi : core) {
    auto f = double(phi.twoj() + 1) * phi.occ_frac;
    for (auto i = 0ul; i < rgrid->num_points; i++) {
      rho[i] += f * (phi.f[i] * phi.f[i] + phi.g[i] * phi.g[i]);
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
    const double rmax, const int stride, const bool include_G,
    const std::vector<double> &lambdas, const std::vector<double> &fk,
    const std::string &fname, const bool FeynmanQ, const bool ScreeningQ,
    const bool holeParticleQ, const int lmax, const bool GreenBasis,
    const bool PolBasis, const double omre, double w0, double wratio) {
  if (valence.empty())
    return;

  // Form list of energies for each kappa:
  std::vector<double> en_list_kappa{};
  if (form_matrix) {
    const auto max_ki =
        std::max_element(cbegin(valence), cend(valence), DiracSpinor::comp_ki)
            ->k_index();
    // for each kappa, find lowest valence (or basis) state energy:
    for (int ki = 0; ki <= max_ki; ki++) {
      // XXX Note: this assumes sorted (for each kappa)! Almost certainly true.
      const auto find_ki = [=](const auto &a) { return a.k_index() == ki; };
      auto vki = std::find_if(cbegin(valence), cend(valence), find_ki);
      // if (vki == cend(valence)) {
      //   vki = std::find_if(cbegin(excited), cend(excited), find_ki);
      // }
      if (vki != cend(valence) /*&& vki != cend(excited)*/) {
        en_list_kappa.push_back(vki->en);
      } else {
        en_list_kappa.push_back(0.0);
      }
    }
  }

  const auto method =
      FeynmanQ ? MBPT::Method::Feynman : MBPT::Method::Goldstone;
  const auto sigp = MBPT::Sigma_params{
      method, nmin_core, include_G, lmax,       GreenBasis,    PolBasis,
      omre,   w0,        wratio,    ScreeningQ, holeParticleQ, fk};
  const auto subgridp = MBPT::rgrid_params{r0, rmax, std::size_t(stride)};

  // Correlaion potential matrix:
  switch (method) {
  case MBPT::Method::Goldstone:
    m_Sigma = std::make_unique<MBPT::GoldstoneSigma2>(
        m_pHF.get(), basis, sigp, subgridp, en_list_kappa, fname);
    break;
  case MBPT::Method::Feynman:
    m_Sigma = std::make_unique<MBPT::FeynmanSigma>(
        m_pHF.get(), basis, sigp, subgridp, en_list_kappa, fname);
    break;
  }

  m_Sigma->scale_Sigma(lambdas);
}

//******************************************************************************
void Wavefunction::hartreeFockBrueckner(const bool print) {
  if (!m_pHF || !m_Sigma) {
    std::cerr << "WARNING 62: Cant call hartreeFockValence before "
                 "hartreeFockCore\n";
    return;
  }
  m_pHF->solveBrueckner(&valence, *(m_Sigma.get()), print);
}

//******************************************************************************
void Wavefunction::fitSigma_hfBrueckner(
    const std::string &valence_list, const std::vector<double> &fit_energies) {
  std::cout << "\nFitting Sigma for lowest valence states:\n";

  std::vector<double> lambdas(fit_energies.size(), 1.0);
  valence.clear();
  hartreeFockValence(valence_list, false);
  const auto hf_val = valence; // copy
  for (int i = 0; true; i++) {
    valence.clear();
    hartreeFockValence(valence_list, false);
    m_Sigma->scale_Sigma(lambdas);
    hartreeFockBrueckner(false);

    auto max_eps = 0.0;
    for (auto ki = 0ul; ki < fit_energies.size(); ++ki) {
      auto match_ki = [=](const auto &v) { return v.k_index() == int(ki); };
      auto e_exp = fit_energies[ki];
      auto l_0 = lambdas[ki];
      // find lowest valence state for this kappa
      auto v_ki = std::find_if(cbegin(valence), cend(valence), match_ki);
      if (v_ki == cend(valence))
        continue;
      const auto vhf_ki = std::find_if(cbegin(hf_val), cend(hf_val), match_ki);
      const auto ehf = vhf_ki->en;
      const auto e0 = v_ki->en;
      // E0 = E_hf + l_0 * sigma
      // Eexp = E_hf + l * sigma
      // => l = l0 * (1 + R)
      // R = (Eexp - E0)/(E0-Ehf)
      const auto r = (e_exp - e0) / (e0 - ehf);
      const auto l = l_0 * (1.0 + r);
      lambdas[ki] = std::clamp(l, 0.5, 1.5);
      const auto eps = std::abs((e_exp - e0) / e_exp);
      max_eps = std::max(max_eps, eps);
    }
    if (max_eps < 1.0e-8 || i == 15) {
      std::cout << "converged to: " << max_eps << " [" << i << "]\n";
      break;
    }
  }
  for (const auto &l : lambdas) {
    std::cout << l << ", ";
  }
  std::cout << "\n";
  valence.clear();
  hartreeFockValence(valence_list, false);
  m_Sigma->scale_Sigma(lambdas);
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
      e0 = (v.en + delta);
    printf("%6s| %9.6f  %+9.6f  %+9.6f | %9.6f = %8.1f  %7.1f\n",
           v.symbol().c_str(), v.en, delta, delta2, (v.en + delta),
           (v.en + delta) * cm, (v.en + delta - e0) * cm);
    if (std::abs(delta / v.en) > 0.2)
      std::cout << "      *** Warning: delta too large?\n";
  }
}

//******************************************************************************
std::vector<double> Wavefunction::get_Vlocal(int l) const {
  return qip::add(vnuc, vdir, vrad.get_Hel(l));
}
//******************************************************************************
const std::vector<double> &Wavefunction::get_Hmag(int l) const {
  return vrad.get_Hmag(l);
}
