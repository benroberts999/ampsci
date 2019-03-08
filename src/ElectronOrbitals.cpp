// class ElectronOrbitals::
#include "ElectronOrbitals.h"
#include "ADAMS_solveLocalBS.h"
#include "ATI_atomInfo.h"
#include "DiracSpinor.h"
#include "FPC_physicalConstants.h"
#include "Grid.h"
#include "Nucleus.h"
#include "NumCalc_quadIntegrate.h"
#include <algorithm> //for sort
#include <cmath>
#include <gsl/gsl_sf_fermi_dirac.h> //For fermiNucleus
#include <sstream>
#include <string>
#include <vector>
/*
 ###= To Do ###=
 * Write out to disk

*/

//******************************************************************************
ElectronOrbitals::ElectronOrbitals(int in_z, int in_a, int in_ngp, double rmin,
                                   double rmax, double var_alpha)
    : rgrid(rmin, rmax, (size_t)in_ngp, GridType::loglinear, 3.5),
      m_alpha(FPC::alpha * var_alpha), m_Z(in_z),
      m_A((in_a < 0) ? ATI::defaultA(m_Z) : in_a)
//
{
  // Use Fermi nucleus by default, unless A=0 is given
  if (m_A > 0)
    formNuclearPotential(NucleusType::Fermi);
  else
    formNuclearPotential(NucleusType::zero);
}

//******************************************************************************
int ElectronOrbitals::solveLocalDirac(int n, int k, double e_a, int log_dele_or,
                                      bool iscore)
/*
Uses ADAMS::solveDBS to solve Dirac Eqn for local potential (Vnuc + Vdir)
*/
{
  orbitals.emplace_back(DiracSpinor{n, k, rgrid.ngp});

  // Fill V(r) with nulcear + DIRECT part of electron potential
  // nb: for exchange part, need to use reSolveDirac()
  std::vector<double> v_a = vnuc;
  if (vdir.size() != 0) {
    for (size_t i = 0; i < rgrid.ngp; i++)
      v_a[i] += vdir[i];
  }

  // Solve local dirac Eq:

  auto &psi = orbitals.back();
  psi.en = e_a == 0 ? enGuessVal(n, k) : e_a;
  ADAMS::solveDBS(psi, v_a, rgrid, m_alpha, log_dele_or);

  auto this_index = stateIndexList.size();
  stateIndexList.push_back(this_index);
  // std::cerr << "this-index = " << this_index << "\n";
  if (iscore)
    coreIndexList.push_back(this_index);
  else
    valenceIndexList.push_back(this_index);

  return 0;
}

//******************************************************************************
int ElectronOrbitals::reSolveDirac(DiracSpinor &psi, double e_a,
                                   const std::vector<double> &vex,
                                   int log_dele_or)
/*
"Re"solves dirac eqaution. Use this to re-solve for same state.
Over-rides existing solution.
If no e_a is given, will use the existing one!
(Usually, a better guess should be given, using P.T.)
Note: optionally takes in exchange potential! (see overloaded above)
Note: Uses the "dodgy" re-scaled exchange potenital:
Vex\psi_a = sum_b vex_a psi_b -> [sum_b vex_a (psi_b/psi_a)] psi_a
so, here, vex = [sum_b vex_a (psi_b/psi_a)]
This is not ideal..
*/
{

  // XXX this is inneficient. Fine, except for HF. THen, slow! ?
  std::vector<double> v_a = vnuc;
  if (vdir.size() != 0)
    for (size_t i = 0; i < rgrid.ngp; i++)
      v_a[i] += vdir[i];
  if (vex.size() != 0)
    for (size_t i = 0; i < rgrid.ngp; i++)
      v_a[i] += vex[i];

  // auto &psi = orbitals[i];
  if (e_a != 0)
    psi.en = e_a;
  ADAMS::solveDBS(psi, v_a, rgrid, m_alpha, log_dele_or);

  return 0;
}

//******************************************************************************
int ElectronOrbitals::reSolveDirac(DiracSpinor &psi, double e_a,
                                   int log_dele_or)
/*Overloaded version; see below
This one doesn't have exchange potential
*/
{
  std::vector<double> dummy_vex;
  return reSolveDirac(psi, e_a, dummy_vex, log_dele_or);
}

//******************************************************************************
double ElectronOrbitals::radialIntegral(const DiracSpinor &psi_a,
                                        const DiracSpinor &psi_b,
                                        Operator op) const {
  std::vector<double> empty_vec;
  return radialIntegral(psi_a, psi_b, empty_vec, op);
}

//******************************************************************************
double ElectronOrbitals::radialIntegral(const DiracSpinor &psi_a,
                                        const DiracSpinor &psi_b,
                                        const std::vector<double> &vint,
                                        Operator op) const
/*
Split into dPsi later??
Better to split, so that you can chain operators!
But: will always be a copy in that case; optimise for unity?
(Just take )
*/
{
  // check that a and b are OK!
  // if (a >= (int)f.size() || b >= (int)f.size())
  // std::cout << "\nFail 97 in EO radInt. Invalid state\n";

  // int pinf = std::min(psi_a.pinf, psi_b.pinf);
  // pinf += (int)NumCalc::Nquad + 10;
  // if (pinf >= rgrid.ngp)
  //   pinf = rgrid.ngp - 1;
  // int pinf = rgrid.ngp; // std::max(psi_a.pinf, psi_b.pinf);
  // int pinf = (int)(0.5 * (psi_a.pinf + psi_b.pinf));
  // int pinf = std::max(psi_a.pinf, psi_b.pinf);
  int pinf = 0; // --> ngp

  std::vector<double> fprime_tmp;
  std::vector<double> gprime_tmp;

  // This is to avoid doing a copy..not sure if it worked, or worth it?
  const std::vector<double> *fprime;
  const std::vector<double> *gprime;

  // Is this a good way? Confusing?
  int ig_sign = 1;

  switch (op) {
  case Operator::unity:
    fprime = &psi_b.f;
    gprime = &psi_b.g;
    break;
  case Operator::gamma0:
    fprime = &psi_b.f;
    gprime = &psi_b.g;
    ig_sign = -1; // XXX CHECK! XXX
    break;
  case Operator::gamma5:
    fprime = &psi_b.g;
    gprime = &psi_b.f;
    ig_sign = -1;
    break;
  case Operator::dr:
    fprime_tmp = NumCalc::derivative(psi_b.f, rgrid.drdu, rgrid.du);
    gprime_tmp = NumCalc::derivative(psi_b.g, rgrid.drdu, rgrid.du);
    fprime = &fprime_tmp;
    gprime = &gprime_tmp;
    break;
  case Operator::dr2:
    fprime_tmp = NumCalc::derivative(psi_b.f, rgrid.drdu, rgrid.du, 2);
    gprime_tmp = NumCalc::derivative(psi_b.g, rgrid.drdu, rgrid.du, 2);
    fprime = &fprime_tmp;
    gprime = &gprime_tmp;
    break;
  default:
    std::cout << "\nError 103 in EO radialInt: unknown operator\n";
    return 0;
  }

  double Rf = 0, Rg = 0;
  if (vint.size() == 0) {
    Rf = NumCalc::integrate(psi_a.f, *fprime, rgrid.drdu, 1, 0, pinf);
    Rg = NumCalc::integrate(psi_a.g, *gprime, rgrid.drdu, 1, 0, pinf);
  } else {
    Rf = NumCalc::integrate(psi_a.f, vint, *fprime, rgrid.drdu, 1, 0, pinf);
    Rg = NumCalc::integrate(psi_a.g, vint, *gprime, rgrid.drdu, 1, 0, pinf);
  }

  return (Rf + ig_sign * Rg) * rgrid.du;
}

//******************************************************************************
int ElectronOrbitals::determineCore(std::string str_core_in)
/*
Takes in a string list for the core configuration, outputs an int list
Takes in previous closed shell (noble), + 'rest' (or just the rest)
E.g:
  Core of Cs: Xe
  Core of Gold: Xe 4f14 5d10
'rest' is in form nLm : n=n, L=l, m=number of electrons in that nl shell.
NOTE: Only works up to n=9, and l=5 [h]
*/
{
  // Move comma-seperated string into an array (vector)
  std::vector<std::string> str_core;
  std::stringstream ss(str_core_in);
  while (ss.good()) {
    std::string substr;
    getline(ss, substr, ',');
    str_core.push_back(substr);
  }

  if (str_core.size() == 0)
    return 0; //?

  std::string NobelGas = str_core[0];
  num_core_shell = ATI::getCoreConfig(NobelGas);

  // Giving a nobel gas is optional:
  int ibeg = (num_core_shell.size() == 0) ? 0 : 1;

  for (size_t i = ibeg; i < str_core.size(); i++) {

    // Parse string, determine config for this term
    int n = std::stoi(str_core[i].substr(0, 1));
    int m = std::stoi(str_core[i].substr(2));
    std::string strl = str_core[i].substr(1, 1);
    int l = ATI::symbol_to_l(strl);

    // Check if this term is valid
    if (m > 4 * l + 2 || l + 1 > n) {
      std::cerr << "FAIL 281 EO: invalid core term: " << str_core[i] << "\n";
      return 2;
    }

    std::vector<int> core_ex;
    // Form int list for this term:
    for (int in = 0; in <= n; in++) {
      for (int il = 0; il < in; il++) {
        if (in == n && il == l)
          core_ex.push_back(m);
        else
          core_ex.push_back(0);
      }
    }

    // Merge this term with the existing core:
    auto size = std::max(core_ex.size(), num_core_shell.size());
    core_ex.resize(size);
    num_core_shell.resize(size);

    for (size_t j = 0; j < num_core_shell.size(); j++) {
      num_core_shell[j] += core_ex[j];
      if (num_core_shell[j] > 4 * ATI::core_l[j] + 2)
        return 2; // check if valid
    }
  }

  // Count number of electrons in the core
  num_core_electrons = 0;
  for (int &num : num_core_shell)
    num_core_electrons += num;

  return 0;
}

//******************************************************************************
bool ElectronOrbitals::isInCore(int n, int k) const
/*
Checks if given state is in the core.
*/
{
  for (auto i : coreIndexList)
    if (n == orbitals[i].n && k == orbitals[i].k)
      return true;
  return false;
}

//******************************************************************************
size_t ElectronOrbitals::getStateIndex(int n, int k) const {
  // auto &state_list = forceVal ? valenceIndexList : stateIndexList;
  for (auto i : stateIndexList)
    if (n == orbitals[i].n && k == orbitals[i].k)
      return i;
  std::cerr << "\nFAIL 290 in EO: Couldn't find state nk=" << n << k << "\n";
  return (int)stateIndexList.size(); // this is an invalid index!
}

//******************************************************************************
int ElectronOrbitals::maxCore_n(int ka) const
/*
Returns the largest n for states with kappa = ka in the core
Note: ka is optional input; if none given, will be 0 (& not used)
(used for energy guesses)
Note: if you give it l instead of kappa, still works!
*/
{
  int max_n = 0;
  for (auto i : coreIndexList) {
    if (orbitals[i].k != ka && ka != 0)
      continue;
    if (orbitals[i].n > max_n)
      max_n = orbitals[i].n;
  }
  return max_n;
}

//******************************************************************************
int ElectronOrbitals::solveInitialCore(std::string str_core, int log_dele_or)
/*
Solves the Dirac eqn for each state in the core
Only for local potential (direct part)
HartreeFockClass.cpp has routines for Hartree Fock
*/
{

  int core_ok = determineCore(str_core);
  if (core_ok == 2) {
    std::cout << "Problem with core: " << str_core << "\n";
    return core_ok;
  }

  for (size_t i = 0; i < num_core_shell.size(); i++) {
    int num = num_core_shell[i];
    if (num == 0)
      continue;
    int n = ATI::core_n[i];
    int l = ATI::core_l[i];
    double en_a = enGuessCore(n, l);
    int k1 = l; // j = l-1/2
    if (k1 != 0) {
      solveLocalDirac(n, k1, en_a, log_dele_or, true);
      en_a = 0.95 * orbitals.back().en;
      if (en_a > 0)
        en_a = enGuessCore(n, l);
    }
    int k2 = -(l + 1); // j=l+1/2
    solveLocalDirac(n, k2, en_a, log_dele_or, true);
  }
  auto num_core_states = orbitals.size(); // store number of states in core

  // occupancy fraction for each core state (avg of Non-rel states!):
  for (size_t i = 0; i < num_core_states; i++) {
    int n = orbitals[i].n;
    int ka = orbitals[i].k;
    int l = ATI::l_k(ka);
    // Find the correct core list index (to determine filling factor):
    auto ic = num_core_shell.size();
    for (size_t j = 0; j < num_core_shell.size(); j++) {
      if (n == ATI::core_n[j] && l == ATI::core_l[j]) {
        ic = j;
        break;
      }
    }
    if (ic == num_core_shell.size()) {
      std::cout << "FAIL 254 in ElectronOrbitals:solveInitialCore\n";
      return 2;
    }
    orbitals[i].occ_frac = double(num_core_shell[ic]) / (4 * l + 2);
  }

  return 0;
}

//******************************************************************************
void ElectronOrbitals::orthonormaliseOrbitals(int num_its)
/*
Forces ALL orbitals to be orthogonal to each other, and normal
Note: workes best if run twice!
|a> ->  |a> - \sum_{b!=a} |b><b|a>
Then:
|a> -> |a> / <a|a>
c_ba = c_ab = <a|b>
num_its is optional parameter. Repeats that many times!
Note: I force all orthog to each other - i.e. double count.
{force <2|1>=0 and then <1|2>=0}
Would be 2x faster not to do this - but that would treat some orbitals special!
Hence factor of 0.5
Note: For HF, should never be called after core is frozen!

XXX Note: This allows wfs to extend past pinf!
==> This causes the possible orthog issues..

*/
{
  size_t Ns = orbitals.size();
  std::vector<std::vector<double>> c_ab(Ns, std::vector<double>(Ns));

  // Calculate c_ab = <a|b>  [only for b>a -- symmetric]
  for (auto a : stateIndexList) {
    for (auto b = a + 1; b < Ns; b++) {
      if (orbitals[a].k != orbitals[b].k)
        continue;
      c_ab[a][b] = 0.5 * radialIntegral(orbitals[a], orbitals[b]);
    }
  }

  // Orthogonalise orbitals:
  for (auto a : stateIndexList) {
    for (auto b : stateIndexList) {
      if (a == b)
        continue;
      if (orbitals[a].k != orbitals[b].k)
        continue;
      // c_ab = c_ba : only calc'd half:
      double cab = (a < b) ? c_ab[a][b] : c_ab[b][a];
      for (size_t ir = 0; ir < rgrid.ngp; ir++) {
        orbitals[a].f[ir] -= cab * orbitals[b].f[ir];
        orbitals[a].g[ir] -= cab * orbitals[b].g[ir];
      }
    }
  }

  for (auto &psi : orbitals) {
    double A = radialIntegral(psi, psi);
    double norm = 1. / sqrt(A);
    for (auto &fa_r : psi.f)
      fa_r *= norm;
    for (auto &ga_r : psi.g)
      ga_r *= norm;
  }

  // If necisary: repeat
  if (num_its > 1)
    orthonormaliseOrbitals(num_its - 1);
}

//******************************************************************************
void ElectronOrbitals::orthonormaliseValence(DiracSpinor &psi_v, int num_its,
                                             bool core_only) const
/*
Force given valence orbital to be orthogonal to
  a) all core orbitals
  b) all _lower_ valence orbitals (lower, by state index)
After the core is 'frozen', don't touch core orbitals!
|v> --> |v> - sum_c |c><c|v>
note: here, c denotes core orbitals + valence orbitals with c<v
if core_only=true, will only orthog phi_v against core orbitals
(only need to do this part before generating exchange potential!)
NB: Is it ok that this is const? Not sure... it's strange way to do it?....
*/
{
  auto num_states_below =
      core_only ? coreIndexList.size() : getStateIndex(psi_v.n, psi_v.k);

  // Calculate the coeficients <c|v> = A_cv
  std::vector<double> A_vc(num_states_below);
  for (size_t ic = 0; ic < (size_t)num_states_below; ic++) {
    if (psi_v.k != orbitals[ic].k)
      continue;
    A_vc[ic] = radialIntegral(psi_v, orbitals[ic]); // no 0.5 here
  }

  // Orthogonalise:
  for (size_t ic = 0; ic < (size_t)num_states_below; ic++) {
    if (psi_v.k != orbitals[ic].k)
      continue;
    double Avc = A_vc[ic];
    for (size_t ir = 0; ir < psi_v.pinf; ir++) {
      // Probably an algorithm for this!
      psi_v.f[ir] -= Avc * orbitals[ic].f[ir];
      psi_v.g[ir] -= Avc * orbitals[ic].g[ir];
    }
  }

  // Re-normalise the valence orbital:
  double A = radialIntegral(psi_v, psi_v);
  double norm = 1. / sqrt(A);
  for (auto &fv_r : psi_v.f)
    fv_r *= norm;
  for (auto &gv_r : psi_v.g)
    gv_r *= norm;

  // If necisary: repeat
  if (num_its > 1)
    orthonormaliseValence(psi_v, num_its - 1, core_only);
}

//******************************************************************************
double ElectronOrbitals::enGuessCore(int n, int l) const
/*
Private
Energy guess for core states. Not perfect, good enough
tot_el = total electrons BELOW
num = num electrons in THIS shell
*/
{

  int tot_el = 0;
  int num = 0;
  for (size_t i = 0; i < num_core_shell.size(); i++) {
    if (l == ATI::core_l[i] && n == ATI::core_n[i]) {
      num = num_core_shell[i];
      break;
    }
    tot_el += num_core_shell[i];
  }

  // effective Z (for energy guess) -- not perfect!
  double Zeff = double(m_Z - tot_el - num);
  if (l == 1)
    Zeff = 1. + double(m_Z - tot_el - 0.5 * num);
  if (l == 2)
    Zeff = 1. + double(m_Z - tot_el - 0.5 * num);
  if (Zeff < 1.)
    Zeff = 1.;

  double en_a = -0.5 * pow(Zeff / n, 2);
  if (n > 1)
    en_a *= 0.5;

  if (n == maxCore_n() - 1)
    en_a *= 1.25;
  if (n == maxCore_n())
    en_a *= 1.5;
  if (n == maxCore_n() + 1 && l == 0 && num == 2)
    en_a = -0.5 * pow(Zeff, 2);

  return en_a;
}

//******************************************************************************
double ElectronOrbitals::enGuessVal(int n, int ka) const
/*Energy guess for valence states. Not perfect, good enough*/
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
      rc = Nucleus::approximate_rc(m_A); // update?
    // XXX rms and half-charge-density not the same! Fix XXX
    vnuc = Nucleus::fermiNuclearPotential(m_Z, t, rc, rgrid.r);
    break;
  case NucleusType::spherical:
    if (rc == 0)
      rc = Nucleus::approximate_rc(m_A); // update?
    vnuc = Nucleus::sphericalNuclearPotential(m_Z, rc, rgrid.r);
    break;
  case NucleusType::zero:
    vnuc = Nucleus::sphericalNuclearPotential(m_Z, 0., rgrid.r);
    break;
  default:
    std::cerr << "\nFail EO:755 - invalid nucleus type?\n";
  }
}

//------------------------------------------------------------------------------
void ElectronOrbitals::sortedEnergyList(std::vector<int> &sort_list,
                                        bool do_sort) const
/*
Outouts a list of integers corresponding to the states
sorted by energy (lowest energy first)
*/
{

  sort_list.clear();

  if (!do_sort) {
    for (size_t i = 0; i < orbitals.size(); i++)
      sort_list.push_back((int)i);
    return;
  }

  std::vector<std::vector<double>> t_en;
  for (size_t i = 0; i < orbitals.size(); i++) {
    t_en.push_back({orbitals[i].en, (double)i + 0.1});
    //+0.1 to prevent rounding error when going from double -> int
  }

  // Function pointer: sort 2D vector by first col!
  bool (*sortCol)(const std::vector<double> &v1,
                  const std::vector<double> &v2) =
      [](const std::vector<double> &v1, const std::vector<double> &v2) {
        return v1[0] > v2[0];
      };

  std::sort(t_en.rbegin(), t_en.rend(), sortCol);

  for (size_t i = 0; i < orbitals.size(); i++)
    sort_list.push_back((int)t_en[i][1]);
}
