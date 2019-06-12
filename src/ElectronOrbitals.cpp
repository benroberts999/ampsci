// class ElectronOrbitals::
#include "ElectronOrbitals.hpp"
#include "ADAMS_bound.hpp"
#include "ATI_atomInfo.hpp"
#include "DiracSpinor.hpp"
#include "FPC_physicalConstants.hpp"
#include "Grid.hpp"
#include "Nucleus.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <algorithm> //for sort
#include <cmath>
#include <gsl/gsl_sf_fermi_dirac.h> //For fermiNucleus
#include <sstream>
#include <string>
#include <vector>

//******************************************************************************
ElectronOrbitals::ElectronOrbitals(int in_z, int in_a, int in_ngp, double rmin,
                                   double rmax, double var_alpha)
    : rgrid(rmin, rmax, (std::size_t)in_ngp, GridType::loglinear, 3.5),
      m_alpha(FPC::alpha * var_alpha), m_Z(in_z),
      m_A((in_a < 0) ? ATI::defaultA(m_Z) : in_a) {
  // Use Fermi nucleus by default, unless A=0 is given
  if (m_A > 0)
    formNuclearPotential(NucleusType::Fermi);
  else
    formNuclearPotential(NucleusType::zero);
}

//******************************************************************************
void ElectronOrbitals::solveDirac(DiracSpinor &psi, double e_a,
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

  // XXX this is inneficient. Fine, except for HF. THen, slow! ?
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

  return;
}
//------------------------------------------------------------------------------
void ElectronOrbitals::solveDirac(DiracSpinor &psi, double e_a,
                                  int log_dele_or) const
// Overloaded version; see above
// This one doesn't have exchange potential
{
  std::vector<double> empty_vec;
  return solveDirac(psi, e_a, empty_vec, log_dele_or);
}

//******************************************************************************
void ElectronOrbitals::determineCore(const std::string &str_core_in)
// Takes in a string list for the core configuration, outputs an int list
// Takes in previous closed shell (noble), + 'rest' (or just the rest)
// E.g:
//   Core of Cs: Xe
//   Core of Gold: Xe 4f14 5d10
// 'rest' is in form nLm : n=n, L=l, m=number of electrons in that nl shell.
// NOTE: Only works up to n=9, and l=5 [h] XXX
{

  // If there's a 'Noble-Gas' term, replace it with full config
  // Otherwise, 'first-term' remains unchanges
  auto found = str_core_in.find(",");
  if (found > str_core_in.length())
    found = str_core_in.length();
  auto first_term = str_core_in.substr(0, found);
  auto rest = str_core_in.substr(found);
  auto str_core = ATI::coreConfig(first_term) + rest;

  // Move comma-seperated string into an array (vector)
  std::vector<std::string> term_str_list;
  {
    std::stringstream ss(str_core);
    while (ss.good()) {
      std::string substr;
      getline(ss, substr, ',');
      term_str_list.push_back(substr);
    }
  }

  bool bad_core = false;
  for (const auto &term : term_str_list) {
    // Parse string, determine config for this term

    bool term_ok = true;
    if (term.size() < 3)
      term_ok = false;

    int n{0}, m{0};
    try {
      // xxx here: nlm: n must be <9. l must be 1 digit. ok?
      n = std::stoi(term.substr(0, 1));
      m = std::stoi(term.substr(2));
    } catch (...) {
      term_ok = false;
    }
    std::string strl = term.substr(1, 1);
    int l = ATI::symbol_to_l(strl);

    // Check if this term is valid
    if (l + 1 > n || l < 0 || n < 1)
      term_ok = false;

    if (!term_ok) {
      std::cout << "Problem with core: " << str_core_in << "\n";
      std::cerr << "invalid core term: " << term << "\n";
      std::abort();
    }

    std::vector<int> single_core_term; //'extra' core parts
    // Form int list for this term:
    // Note has entry for every term (most of them zero)
    // Ineficient, but doesn't matter
    for (int in = 0; in <= n; in++) {
      for (int il = 0; il < in; il++) {
        if (in == n && il == l)
          single_core_term.push_back(m);
        else
          single_core_term.push_back(0);
      }
    }

    // Merge this term with the existing core:
    auto size = std::max(single_core_term.size(), num_core_shell.size());
    single_core_term.resize(size);
    num_core_shell.resize(size);
    for (std::size_t j = 0; j < num_core_shell.size(); j++) {
      num_core_shell[j] += single_core_term[j];
      if (num_core_shell[j] > 4 * ATI::core_l[j] + 2 || num_core_shell[j] < 0)
        bad_core = true;
    }
  }

  if (bad_core) {
    std::cout << "Problem with core: " << str_core_in << " = \n";
    for (std::size_t j = 0; j < num_core_shell.size(); j++) {
      auto num = num_core_shell[j];
      auto n = ATI::core_n[j];
      auto l = ATI::core_l[j];
      if (num == 0)
        continue;
      if (num > 4 * l + 2 || l + 1 > n || num < 0)
        std::cout << " **";
      std::cout << n << ATI::l_symbol(l) << num << ",";
    }
    std::cout << "\n";
    std::cout << "In this house, we obey the Pauli exclusion principle!\n";
    std::abort();
  }

  // Count number of electrons in the core
  num_core_electrons = 0;
  for (int &num : num_core_shell)
    num_core_electrons += num;

  if (num_core_electrons > m_Z) {
    std::cout << "Problem with core: " << str_core_in << "\n";
    std::cout << "Too many electrons: N_core=" << num_core_electrons
              << ", Z=" << m_Z << "\n";
    std::abort();
  }

  // store full core config list
  m_core_string = str_core;
  // XXX Would be nice to "merge" any duplicates (but retain input order)
  // i.e. so that "...,5p6,...,5p-1" --> "...,5p5,..."

  return;
}

//******************************************************************************
bool ElectronOrbitals::isInCore(int n, int k) const
// Checks if given state is in the core.
{
  for (auto &phi : core_orbitals) {
    if (n == phi.n && k == phi.k)
      return true;
  }
  return false;
}

//******************************************************************************
std::size_t ElectronOrbitals::getStateIndex(int n, int k) const {
  for (auto &tmp_orbitals : {core_orbitals, valence_orbitals}) {
    // XXX How the hell does this work?
    // Is is making a copy of core_orbitals + valence_orbitals ?
    // If so, have to be careful not to use this deep in a loop!!!
    for (std::size_t i = 0; i < tmp_orbitals.size(); i++) {
      auto &phi = tmp_orbitals[i];
      if (n == phi.n && k == phi.k)
        return i;
    }
  }
  std::cerr << "\nFAIL 290 in EO: Couldn't find state n,k=" << n << "," << k
            << "\n";
  std::abort();
}

//******************************************************************************
int ElectronOrbitals::maxCore_n(int ka) const
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
int ElectronOrbitals::solveInitialCore(std::string str_core, int log_dele_or)
// Solves the Dirac eqn for each state in the core
// Only for local potential (direct part)
// HartreeFockClass.cpp has routines for Hartree Fock
{
  if (core_orbitals.size() > 0) {
    std::cerr << "Fail 254 in ElectronOrbitals:solveInitialCore: States "
                 "already exist! "
              << core_orbitals.size() << "\n";
    std::abort();
  }

  determineCore(str_core);

  for (std::size_t i = 0; i < num_core_shell.size(); i++) {
    int num = num_core_shell[i];
    if (num == 0)
      continue;
    int n = ATI::core_n[i];
    int l = ATI::core_l[i];
    double en_a = enGuessCore(n, l);
    int k1 = l; // j = l-1/2
    if (k1 != 0) {
      core_orbitals.emplace_back(DiracSpinor{n, k1, rgrid});
      solveDirac(core_orbitals.back(), en_a, log_dele_or);
      en_a = 0.95 * core_orbitals.back().en;
      if (en_a > 0)
        en_a = enGuessCore(n, l);
    }
    int k2 = -(l + 1); // j=l+1/2
    core_orbitals.emplace_back(DiracSpinor{n, k2, rgrid});
    solveDirac(core_orbitals.back(), en_a, log_dele_or);
  }

  // occupancy fraction for each core state (avg of Non-rel states!):
  for (auto &phi : core_orbitals) {
    int n = phi.n;
    int l = phi.l();
    // Find the correct core list index (to determine filling factor):
    auto ic = num_core_shell.size();
    for (std::size_t j = 0; j < num_core_shell.size(); j++) {
      if (n == ATI::core_n[j] && l == ATI::core_l[j]) {
        ic = j;
        break;
      }
    }
    if (ic == num_core_shell.size()) {
      std::cout << "FAIL 254 in ElectronOrbitals:solveInitialCore\n";
      return 2;
    }
    phi.occ_frac = double(num_core_shell[ic]) / (4 * l + 2);
  }

  return 0;
}

//******************************************************************************
void ElectronOrbitals::solveNewValence(int n, int k, double en_a,
                                       int log_dele_or)
// Update to take a list ok nken's ?
{

  valence_orbitals.emplace_back(DiracSpinor{n, k, rgrid});

  // XXX Do this again? Or just pass to solveDirac?
  // Fill V(r) with nulcear + DIRECT part of electron potential
  std::vector<double> v_a = vnuc;
  if (vdir.size() != 0) {
    for (auto i = 0ul; i < rgrid.ngp; i++)
      v_a[i] += vdir[i];
  }

  // Solve local dirac Eq:
  auto &psi = valence_orbitals.back();
  psi.en = en_a == 0 ? enGuessVal(n, k) : en_a;
  ADAMS::solveDBS(psi, v_a, rgrid, m_alpha, log_dele_or);
}

//******************************************************************************
void ElectronOrbitals::orthonormaliseOrbitals(
    std::vector<DiracSpinor> &tmp_orbs, int num_its)
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
  std::size_t Ns = tmp_orbs.size();
  std::vector<std::vector<double>> c_ab(Ns, std::vector<double>(Ns));

  // Calculate c_ab = <a|b>  [only for b>a -- symmetric]
  // for (auto a : tmp_orbs) {
  for (std::size_t a = 0; a < Ns; a++) {
    for (auto b = a + 1; b < Ns; b++) {
      if (tmp_orbs[a].k != tmp_orbs[b].k)
        continue;
      c_ab[a][b] = 0.5 * (tmp_orbs[a] * tmp_orbs[b]);
    }
  }

  // Orthogonalise orbitals:
  for (std::size_t a = 0; a < Ns; a++) {
    for (std::size_t b = 0; b < Ns; b++) {
      if (tmp_orbs[a].k != tmp_orbs[b].k)
        continue;
      if (a == b)
        continue;
      // c_ab = c_ba : only calc'd half:
      double cab = (a < b) ? c_ab[a][b] : c_ab[b][a];
      for (std::size_t ir = 0; ir < tmp_orbs[a].p_rgrid->ngp; ir++) {
        tmp_orbs[a].f[ir] -= cab * tmp_orbs[b].f[ir];
        tmp_orbs[a].g[ir] -= cab * tmp_orbs[b].g[ir];
      }
    }
  }

  for (auto &psi : tmp_orbs) {
    double norm = 1. / sqrt(psi * psi);
    for (auto &fa_r : psi.f)
      fa_r *= norm;
    for (auto &ga_r : psi.g)
      ga_r *= norm;
  }

  // If necisary: repeat
  if (num_its > 1)
    orthonormaliseOrbitals(tmp_orbs, num_its - 1);
}

//******************************************************************************
void ElectronOrbitals::orthonormaliseWrtCore(DiracSpinor &psi_v) const
// Force given orbital to be orthogonal to all core orbitals
// [After the core is 'frozen', don't touch core orbitals!]
// |v> --> |v> - sum_c |c><c|v>
// note: here, c denotes core orbitals
{

  auto Nc = core_orbitals.size();

  // Calculate the coeficients <c|v> = A_cv
  std::vector<double> A_vc; //(num_states_below);
  A_vc.reserve(Nc);
  for (const auto &phi_c : core_orbitals) {
    if (psi_v.k != phi_c.k) {
      A_vc.push_back(0.);
      continue;
    }
    A_vc.emplace_back((psi_v * phi_c)); // no 0.5 here
  }

  // Orthogonalise:
  for (std::size_t ic = 0; ic < Nc; ic++) {
    const auto &psi_c = core_orbitals[ic];
    if (psi_v.k != psi_c.k)
      continue;
    const double Avc = A_vc[ic];
    for (std::size_t ir = 0; ir < psi_v.pinf; ir++) {
      // Probably an algorithm for this!
      psi_v.f[ir] -= Avc * psi_c.f[ir];
      psi_v.g[ir] -= Avc * psi_c.g[ir];
    }
  }

  // Re-normalise the valence orbital:
  const double norm = 1. / sqrt((psi_v * psi_v));
  for (auto &fv_r : psi_v.f)
    fv_r *= norm;
  for (auto &gv_r : psi_v.g)
    gv_r *= norm;
}

//******************************************************************************
double ElectronOrbitals::enGuessCore(int n, int l) const
// Private
// Energy guess for core states. Not perfect, good enough
// tot_el = total electrons BELOW
// num = num electrons in THIS shell
{
  int tot_el = 0;
  int num = 0;
  for (std::size_t i = 0; i < num_core_shell.size(); i++) {
    if (l == ATI::core_l[i] && n == ATI::core_n[i]) {
      num = num_core_shell[i];
      break;
    }
    tot_el += num_core_shell[i];
  }

  // effective Z (for energy guess) -- not perfect!
  double Zeff = (m_Z - tot_el - num);
  if (l == 1) {
    Zeff = 1. + (m_Z - tot_el - 0.5 * num);
  } else if (l == 2) {
    Zeff = 1. + (m_Z - tot_el - 0.5 * num);
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
  } else if (n == maxCore_n() + 1 && l == 0 && num == 2) {
    en_a = -0.5 * pow(Zeff, 2);
  }

  return en_a;
}

//******************************************************************************
double ElectronOrbitals::enGuessVal(int n, int ka) const
// Energy guess for valence states. Not perfect, good enough
{
  int maxn = maxCore_n();
  int l = ATI::l_k(ka);
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
void ElectronOrbitals::formNuclearPotential(NucleusType nucleus_type, double rc,
                                            double t) {
  vnuc.clear();
  switch (nucleus_type) {
  case NucleusType::Fermi:
    if (t == 0)
      t = Nucleus::approximate_t_skin(m_A);
    if (rc == 0)
      rc = Nucleus::approximate_c_hdr(m_A);
    vnuc = Nucleus::fermiNuclearPotential(m_Z, t, rc, rgrid.r);
    break;
  case NucleusType::spherical:
    if (rc == 0) {
      // note: still called rc, but is r_N here!
      rc = Nucleus::approximate_r_rms(m_A);
    }
    t = 0;
    vnuc = Nucleus::sphericalNuclearPotential(m_Z, rc, rgrid.r);
    break;
  case NucleusType::zero:
    rc = 0;
    t = 0;
    vnuc = Nucleus::sphericalNuclearPotential(m_Z, 0., rgrid.r);
    break;
  default:
    std::cerr << "\nFail EO:755 - invalid nucleus type?\n";
  }
  m_c = rc;
  m_t = t;
}

//******************************************************************************
void ElectronOrbitals::printNuclearParams() {
  if (m_c == 0 && m_t == 0) {
    std::cout << "Zero-size nucleus\n";
  } else if (m_t == 0) {
    std::cout << "Spherical nucleus; r_rms = " << m_c << "\n";
  } else {
    std::cout << "Fermi nucleus; r_rms = "
              << Nucleus::rrms_formula_c_t(m_c, m_t) << ", c=" << m_c
              << ", t=" << m_t << "\n";
  }
}

//******************************************************************************
std::vector<std::size_t>
ElectronOrbitals::sortedEnergyList(const std::vector<DiracSpinor> &tmp_orbs,
                                   bool do_sort) const
// Outouts a list of integers corresponding to the states
// sorted by energy (lowest energy first)
{
  std::vector<std::vector<double>> t_en;
  for (std::size_t i = 0; i < tmp_orbs.size(); i++) {
    t_en.push_back({tmp_orbs[i].en, (double)i + 0.1});
    //+0.1 to prevent rounding error when going from double -> int
  }

  // Function pointer: sort 2D vector by first col
  auto sortCol = [](const std::vector<double> &v1,
                    const std::vector<double> &v2) { return v1[0] > v2[0]; };

  if (do_sort)
    std::sort(t_en.rbegin(), t_en.rend(), sortCol);

  // overwrite list with sorted list
  std::vector<std::size_t> sorted_list;
  for (const auto &el : t_en) {
    sorted_list.push_back(std::size_t(el[1]));
  }

  return sorted_list;
}

//******************************************************************************
void ElectronOrbitals::printCore(bool sorted)
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
           phi.en *FPC::Hartree_invcm);
    if (phi.occ_frac < 1.0) {
      printf("     (%4.2f)\n", phi.occ_frac);
    } else {
      std::cout << "\n";
    }
  }
}

//******************************************************************************
void ElectronOrbitals::printValence(bool sorted)
// prints valence orbitals
{
  if (valence_orbitals.size() == 0)
    return;

  std::cout << "Val: state   "
            << "k   Rinf its   eps       En (au)      En (/cm)   En (/cm)\n";

  // Find lowest valence energy:
  double e0 = 0;
  for (auto &phi : valence_orbitals) {
    if (phi.en < e0)
      e0 = phi.en;
  }

  auto index_list = sortedEnergyList(valence_orbitals, sorted);
  for (auto i : index_list) {
    auto &phi = valence_orbitals[i];
    double r_inf = rinf(phi);
    printf("%2i) %7s %2i  %5.1f %2i  %5.0e %13.7f %13.1f", int(i),
           phi.symbol().c_str(), phi.k, r_inf, phi.its, phi.eps, phi.en,
           phi.en *FPC::Hartree_invcm);
    printf(" %10.2f\n", (phi.en - e0) * FPC::Hartree_invcm);
  }
}

//******************************************************************************
std::vector<std::vector<int>>
ElectronOrbitals::listOfStates_nk(int num_val, int la, int lb, bool skip_core)
// Creates a list of states (usually valence states to solve for)
// In form {{n,ka},...}
// Outputs n and l, from min-max l (or from 0-> max l)
// Calulates num_val different n states (for each l)
// This will be number of states above the core (skips states which are in the
// core)
{
  std::vector<std::vector<int>> lst;
  auto l_min = la;
  auto l_max = lb;
  if (lb == 0) {
    l_min = 0;
    l_max = la;
  }

  // XXX a) Add energy guess for valence/core!
  // XXX b) Change to give list of nken's !
  // XXX Then: send this list to solvers

  auto min_ik = ATI::indexFromKappa(-l_min - 1);
  auto max_ik = ATI::indexFromKappa(-l_max - 1);
  for (int ik = min_ik; ik <= max_ik; ik++) {
    auto k = ATI::kappaFromIndex(ik);
    auto l = ATI::l_k(k);
    auto n_min = l + 1;
    for (int n = n_min, count = 0; count < num_val; n++) {
      if (isInCore(n, k) && skip_core)
        continue;
      lst.push_back({n, k});
      ++count;
    }
  }
  return lst;
}
