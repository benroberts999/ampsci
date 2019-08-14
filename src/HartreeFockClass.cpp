#include "HartreeFockClass.hpp"
#include "Adams/Adams_bound.hpp"
#include "CoulombIntegrals.hpp"
#include "DiracOperator.hpp"
#include "DiracSpinor.hpp"
#include "Grid.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include "Operators.hpp"
#include "Parametric_potentials.hpp"
#include "Physics/AtomInfo.hpp"
#include "Physics/Wigner_369j.hpp"
#include "Wavefunction.hpp"
#include <cmath>
#include <vector>
/*
Calculates self-consistent Hartree-Fock potential, including exchange.
Solves all core and valence states.
*/

#define DO_DEBUG false
#if DO_DEBUG
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif // DEBUG

//******************************************************************************
HartreeFock::HartreeFock(HFMethod method, Wavefunction &wf,
                         const std::string &in_core, double eps_HF, double h_d,
                         double g_t)
    //
    : p_wf(&wf), p_rgrid(&wf.rgrid),
      m_cint(Coulomb(wf.rgrid, wf.core_orbitals, wf.valence_orbitals)),
      m_eps_HF([=]() { // can give as log..
        return (std::fabs(eps_HF) < 1) ? eps_HF : std::pow(10, -1 * eps_HF);
      }()),
      m_excludeExchange([=]() {
        return (method == HFMethod::HartreeFock || method == HFMethod::ApproxHF)
                   ? false
                   : true;
      }()), //
      m_method(method)
//
{
  // XXX Update this so that you can create HF class, THEN solve for core later
  // if(wf.Z)
  if (method == HFMethod::Hartree || method == HFMethod::HartreeFock ||
      method == HFMethod::ApproxHF) {
    if (wf.core_orbitals.empty())
      starting_approx_core(in_core, 5);
    m_cint.initialise_core_core();
    appr_vex_core.resize(p_wf->core_orbitals.size(),
                         std::vector<double>(p_rgrid->ngp));
    hartree_fock_core();
  } else {
    starting_approx_core(in_core, 15, method, h_d, g_t);
    m_cint.initialise_core_core();
    appr_vex_core.resize(p_wf->core_orbitals.size(),
                         std::vector<double>(p_rgrid->ngp));
    m_cint.form_core_core();
  }
}

//------------------------------------------------------------------------------
// Overload (to allow new HF object for seperate basis..)
// very hacky temp solution!
HartreeFock::HartreeFock(Wavefunction &wf,
                         const std::vector<DiracSpinor> &val_orbitals,
                         double eps_HF, bool in_ExcludeExchange)
    : p_wf(&wf), p_rgrid(&wf.rgrid),
      m_cint(Coulomb(wf.rgrid, wf.core_orbitals, val_orbitals)),
      m_eps_HF([=]() { // can give as log..
        return (std::fabs(eps_HF) < 1) ? eps_HF : std::pow(10, -1 * eps_HF);
      }()),
      m_excludeExchange(in_ExcludeExchange), //
      m_method(HFMethod::HartreeFock)
// Core must already exist to use this one!
// Call it something else??
{
  m_cint.initialise_core_core();
}

//******************************************************************************
HFMethod HartreeFock::parseMethod(std::string in_method) {
  if (in_method == "HartreeFock")
    return HFMethod::HartreeFock;
  if (in_method == "ApproxHF")
    return HFMethod::ApproxHF;
  if (in_method == "Hartree")
    return HFMethod::Hartree;
  if (in_method == "GreenPRM")
    return HFMethod::GreenPRM;
  if (in_method == "TietzPRM")
    return HFMethod::TietzPRM;
  std::cout << "Warning: HF Method: " << in_method << " ?? Defaulting to HF\n";
  return HFMethod::HartreeFock;
}

//******************************************************************************
void HartreeFock::hartree_fock_core() {

  if (p_wf->core_orbitals.empty()) {
    // If H-like, kill "initial" vdir (Green potential)
    p_wf->vdir = std::vector<double>(p_wf->rgrid.ngp, 0);
    return;
  }

  double eps_target_HF = (m_method == HFMethod::ApproxHF) ? m_eps_HF : 1.0e-5;
  bool do_refine = true;
  if (p_wf->core_orbitals.size() == 1) {
    // in this case, "Approximate" method is "exact" (And green fails???)
    eps_target_HF = m_eps_HF;
    do_refine = false;
  }

  static const double eta1 = 0.35;
  static const double eta2 = 0.7; // this value after 4 its
  // don't include all pts in PT for new e guess:
  static const std::size_t de_stride = 5;

  // initialise 'old' potentials
  auto vdir_old = p_wf->vdir;
  auto vex_old = appr_vex_core;

  // Start the HF itterative procedure:
  int hits = 1;
  double t_eps = 1.0;
  auto t_eps_prev = 1.0;
  double eta = 1.0;
  for (; hits < MAX_HART_ITS; hits++) {
    DEBUG(std::cerr << "HF core it: " << hits << "\n";)
    if (hits == 2)
      eta = eta1;
    else if (hits == 4)
      eta = eta2;
    else if (hits == 16)
      eta = 0.5 * (eta1 + eta2);
    else if (hits == 32)
      eta = eta1;

    // Store old vdir/vex
    vdir_old = p_wf->vdir;
    vex_old = appr_vex_core;

    // Form new v_dir and v_ex:
    m_cint.form_core_core();
    form_vdir(p_wf->vdir, false);
    form_approx_vex_core(appr_vex_core);
    if (hits == 1)
      vex_old = appr_vex_core; // We didn't have old vex before

    for (std::size_t j = 0; j < p_rgrid->ngp; j++) {
      p_wf->vdir[j] = eta * p_wf->vdir[j] + (1. - eta) * vdir_old[j];
      for (std::size_t i = 0; i < p_wf->core_orbitals.size(); i++) {
        appr_vex_core[i][j] =
            eta * appr_vex_core[i][j] + (1. - eta) * vex_old[i][j];
      }
    }

    // Solve Dirac Eq. for each state in core, using Vdir+Vex:
    t_eps = 0;
    for (std::size_t i = 0; i < p_wf->core_orbitals.size(); i++) {
      auto &phi = p_wf->core_orbitals[i];
      double en_old = phi.en;
      // calculate de from PT
      double del_e = 0;
      for (std::size_t j = 0; j < phi.pinf; j += de_stride) {
        double dv = (p_wf->vdir[j] - vdir_old[j]) +
                    (appr_vex_core[i][j] - vex_old[i][j]);
        del_e += dv * phi.f[j] * phi.f[j] * p_rgrid->drdu[j];
      }
      del_e *= p_rgrid->du * de_stride;
      double en_guess = (en_old < -del_e) ? en_old + del_e : en_old;
      p_wf->solveDirac(phi, en_guess, appr_vex_core[i], 6);
      double state_eps = fabs((phi.en - en_old) / en_old);
      // convergance based on worst orbital:
      DEBUG(printf(" --- %2i,%2i: en=%11.5f  HFeps = %.0e;  Adams = %.0e[%2i]  "
                   "(%4i)\n",
                   phi.n, phi.k, phi.en, state_eps, phi.eps, phi.its,
                   (int)phi.pinf);)
      t_eps = (state_eps > t_eps) ? state_eps : t_eps;
    } // core states
    DEBUG(std::cerr << "HF core it: " << hits << ": eps=" << t_eps << "\n\n";
          std::cin.get();)

    // Force all core orbitals to be orthogonal to each other
    p_wf->orthonormaliseOrbitals(p_wf->core_orbitals, 1);
    auto getting_worse = (hits > 20 && t_eps > t_eps_prev);
    auto converged = (t_eps < eps_target_HF);
    if (converged || getting_worse)
      break;
    t_eps_prev = t_eps;
  } // hits
  if (verbose)
    printf("HF core      it:%3i eps=%6.1e              \n", hits, t_eps);

  // Now, re-solve core orbitals with higher precission
  for (std::size_t i = 0; i < p_wf->core_orbitals.size(); i++) {
    p_wf->solveDirac(p_wf->core_orbitals[i], p_wf->core_orbitals[i].en,
                     appr_vex_core[i], 15);
  }
  p_wf->orthonormaliseOrbitals(p_wf->core_orbitals, 2);

  if (!m_excludeExchange && do_refine)
    refine_core_orbitals_exchange();
}

//******************************************************************************
void HartreeFock::solveNewValence(int n, int kappa) {

  p_wf->valence_orbitals.emplace_back(DiracSpinor{n, kappa, p_wf->rgrid});
  // Solve local dirac Eq:
  auto &phi = p_wf->valence_orbitals.back();
  appr_vex_val.emplace_back(std::vector<double>{});
  auto &vexa = appr_vex_val.back();
  solveValence(phi, vexa);
}

//******************************************************************************
void HartreeFock::solveValence(DiracSpinor &phi, std::vector<double> &vexa)
// Solves HF for given orbital phi, in frozen core.
// Does not store vex (must be done outside)
// Can be used to generate a set of virtual/basis orbitals
{
  auto kappa = phi.k;
  int twoJplus1 = AtomInfo::twoj_k(kappa) + 1;
  phi.occ_frac = 1. / twoJplus1;

  double eps_target_HF = (m_method == HFMethod::ApproxHF) ? m_eps_HF : 1.0e-5;

  static const double eta1 = 0.35;
  static const double eta2 = 0.7; // this value after 4 its
  // don't include all pts in PT for new e guess
  static const std::size_t de_stride = 5;

  vexa.clear();
  vexa.resize(p_rgrid->ngp, 0);

  auto vexa_old = vexa;

  int hits = 1;
  double eps = -1, eps_prev = -1;
  double eta = eta1;
  for (; hits < MAX_HART_ITS; hits++) {
    if (hits == 4)
      eta = eta2;

    double en_old = phi.en;
    vexa_old = vexa;

    m_cint.form_core_valence(phi);
    form_approx_vex_a(phi, vexa);

    for (std::size_t i = 0; i < p_rgrid->ngp; i++) {
      vexa[i] = eta * vexa[i] + (1. - eta) * vexa_old[i];
    }
    // Use P.T. to calculate energy change:
    double en_new_guess = 0;
    for (std::size_t i = 0; i < phi.pinf; i += de_stride) {
      en_new_guess +=
          (vexa[i] - vexa_old[i]) * phi.f[i] * phi.f[i] * p_rgrid->drdu[i];
    }
    en_new_guess = en_old + en_new_guess * p_rgrid->du * de_stride;
    // Solve Dirac using new potential:
    p_wf->solveDirac(phi, en_new_guess, vexa, 15);
    eps = fabs((phi.en - en_old) / en_old);
    // Force valence states to be orthogonal to core:
    p_wf->orthonormaliseWrtCore(phi);

    auto getting_worse = (hits > 20 && eps >= eps_prev);
    auto converged = (eps <= eps_target_HF);
    if (converged || getting_worse)
      break;
    eps_prev = eps;
  }
  if (verbose)
    printf("HF val: %2i %2i | %3i eps=%6.1e  en=%11.8f\n", phi.n, kappa, hits,
           eps, phi.en);

  // Re-solve w/ higher precission
  p_wf->solveDirac(phi, phi.en, vexa, 15);
  p_wf->orthonormaliseWrtCore(phi);

  if (!m_excludeExchange && p_wf->core_orbitals.size() != 0)
    refine_valence_orbital_exchange(phi);
}

//******************************************************************************
double HartreeFock::calculateCoreEnergy() const
// Calculates the total HF core energy:
//   E = \sum_a [ja]e_a - 0.5 \sum_(ab) (R^0_abab - \sum_k L^k_ab R^k_abba)
// where:
//   R^k_abcd = Integral [f_a*f_c + g_a*g_c] * v^k_bd
//   R^0_abab is not absymmetric
//   R^k_abba _is_ ab symmetric
{
  double Etot = 0;
  for (std::size_t a = 0; a < p_wf->core_orbitals.size(); a++) {
    const auto &phi_a = p_wf->core_orbitals[a];
    auto tja = phi_a.twoj();

    double E1 = 0, E2 = 0, E3 = 0;
    double xtjap1 = (tja + 1) * phi_a.occ_frac;
    E1 += xtjap1 * phi_a.en;
    for (std::size_t b = 0; b < p_wf->core_orbitals.size(); b++) {
      const auto &phi_b = p_wf->core_orbitals[b];
      auto tjb = phi_b.twoj();
      double xtjbp1 = (tjb + 1) * phi_b.occ_frac;
      auto irmax = std::min(phi_a.pinf, phi_b.pinf);
      auto &v0bb = m_cint.get_y_ijk(phi_b, phi_b, 0);
      double R0f2 = NumCalc::integrate(
          {&phi_a.f, &phi_a.f, &v0bb, &p_rgrid->drdu}, 1.0, 0, irmax);
      double R0g2 = NumCalc::integrate(
          {&phi_a.g, &phi_a.g, &v0bb, &p_rgrid->drdu}, 1.0, 0, irmax);
      E2 += xtjap1 * xtjbp1 * (R0f2 + R0g2);
      // take advantage of symmetry for third term:
      if (b > a)
        continue;
      double y = (a == b) ? 1 : 2;
      int kmin = std::abs(tja - tjb) / 2;
      int kmax = (tja + tjb) / 2;
      auto &vabk = m_cint.get_y_ijk(phi_a, phi_b);
      const auto &L_abk =
          m_cint.get_angular_L_kiakib_k(phi_a.k_index(), phi_b.k_index());
      for (int k = kmin; k <= kmax; k++) {

        if (L_abk[k - kmin] == 0)
          continue;
        int ik = k - kmin;
        double R0f3 =
            NumCalc::integrate({&phi_a.f, &phi_b.f, &vabk[ik], &p_rgrid->drdu});
        double R0g3 =
            NumCalc::integrate({&phi_a.g, &phi_b.g, &vabk[ik], &p_rgrid->drdu});
        E3 += y * xtjap1 * xtjbp1 * L_abk[k - kmin] * (R0f3 + R0g3);
      }
    }
    {
      Etot += E1 - 0.5 * (E2 - E3) * p_rgrid->du; // update running total
    }
  }
  return Etot;
}

//******************************************************************************
void HartreeFock::starting_approx_core(const std::string &in_core,
                                       int log_converge, HFMethod method,
                                       double h_g, double d_t)
// Starting approx for HF. Uses Green parametric
// Later, can put other options if you want.
{

  if (method == HFMethod::GreenPRM) {
    p_wf->vdir = Parametric::GreenPotential(p_wf->Znuc(), p_rgrid->r, h_g, d_t);
  } else if (method == HFMethod::TietzPRM) {
    p_wf->vdir = Parametric::TietzPotential(p_wf->Znuc(), p_rgrid->r, h_g, d_t);
  } else {
    std::cerr << "FAIL 321: Wrong core starting approx method \n";
    std::abort();
  }
  p_wf->solveInitialCore(in_core, log_converge);
}

//******************************************************************************
void HartreeFock::form_vdir(std::vector<double> &vdir, bool re_scale) const
// Forms the direct part of the potential.
// Must call either form_vbb0 or form_vabk_core first!
// Doesn't calculate, assumes m_arr_v_abk_r array exists + is up-to-date
// If re_scale==true, will scale by (N-1)/N. This then given the averaged
// Hartree potential (local, same each state, no exchange).
// re_scale=false by default
{
  for (auto &v_dir : vdir) {
    v_dir = 0;
  }
  double sf = re_scale ? (1. - 1. / p_wf->Ncore()) : 1;
  for (const auto &phi_b : p_wf->core_orbitals) {
    double f = (phi_b.twoj() + 1) * phi_b.occ_frac;
    const auto &v0bb = m_cint.get_y_ijk(phi_b, phi_b, 0);
    for (std::size_t i = 0; i < p_rgrid->ngp; i++) {
      vdir[i] += f * v0bb[i] * sf;
    }
  }
}

//******************************************************************************
void HartreeFock::form_approx_vex_core(
    std::vector<std::vector<double>> &vex) const
// Forms the 2D "approximate" exchange potential for each core state, a.
// NOTE: Must call form_vabk_core first!
// Doesn't calculate, assumes m_arr_v_abk_r array exists + is up-to-date
{
#pragma omp parallel for
  for (std::size_t a = 0; a < p_wf->core_orbitals.size(); a++) {
    form_approx_vex_a(p_wf->core_orbitals[a], vex[a]);
  }
}

//******************************************************************************
void HartreeFock::form_approx_vex_a(const DiracSpinor &phi_a,
                                    std::vector<double> &vex_a) const
// Forms the 2D "approximate" exchange potential for given core state, a.
// Does the a=b case seperately, since it's a little simpler
// Approximate:
// In order to approximate solution to HF equations, I form "local" ex.
// potential
//   [v_ex*psi_a](r) = \sum_b v_ex^(a,b)(r) * psi_b(r)
// v_ex is non-local; cannot write: [v_ex*psi_a](r) =/= v_ex(r)*psi_a(r)
// Instead: define local approx: vex_a
//   vex_a = [v_ex*psi_a](r) *(psi_a/psi_a^2)
//         = \sum_b v_ex^(a,b)(r)*psi_b(r) * (psi_a/psi_a^2)
//         = \sum_b v_ex^(a,b)(r)*(psi_b(r)*psi_a) / psi_a^2
// This vex_a is then a local potential (different for each state!) that can
// be used as an addition to local direct potential to solve Dirac Eq. as
// normal. In theory, this is exact. Clearly, however, there is an issue when
// psi_a is small. Luckily, however, we don't care as much when psi_a is small!
// Also, since v_ex is already small (compared to vdir), we can make good
// approximation. Therefore, I only calculate vex_a when a=b, or when |psi_a|
// > 1.e3 Further, largest part of v_ex is when a=b. In this case, the factor=1
// is exact!
{
  for (auto &va : vex_a) {
    va = 0;
  }

  auto ki_a = phi_a.k_index();
  auto twoj_a = phi_a.twoj();

  bool a_in_coreQ = false;

  if (!m_excludeExchange) {
    for (const auto &phi_b : p_wf->core_orbitals) { // b!=a
      if (phi_b == phi_a) {
        a_in_coreQ = true;
        continue;
      }
      auto tjb = phi_b.twoj();
      double x_tjbp1 = (tjb + 1) * phi_b.occ_frac;
      auto irmax = std::min(phi_a.pinf, phi_b.pinf);
      int kmin = std::abs(twoj_a - tjb) / 2;
      int kmax = (twoj_a + tjb) / 2;
      const auto &vabk = m_cint.get_y_ijk(phi_b, phi_a);

      // hold "fraction" psi_a*psi_b/(psi_a^2):
      std::vector<double> v_Fab(p_rgrid->ngp);
      for (std::size_t i = 0; i < irmax; i++) {
        // This is the approximte part! Divides by psi_a
        if (std::fabs(phi_a.f[i]) < 1.e-3)
          continue;
        double fac_top = phi_a.f[i] * phi_b.f[i] + phi_a.g[i] * phi_b.g[i];
        double fac_bot = phi_a.f[i] * phi_a.f[i] + phi_a.g[i] * phi_a.g[i];
        v_Fab[i] = -1. * x_tjbp1 * fac_top / fac_bot;
      } // r
      const auto &L_ab_k = m_cint.get_angular_L_kiakib_k(ki_a, phi_b.k_index());
      for (int k = kmin; k <= kmax; k++) {
        if (L_ab_k[k - kmin] == 0)
          continue;
        for (std::size_t i = 0; i < irmax; i++) {
          if (v_Fab[i] == 0)
            continue;
          vex_a[i] += L_ab_k[k - kmin] * vabk[k - kmin][i] * v_Fab[i];
        } // r
      }   // k
    }     // b
  }

  // now, do a=b, ONLY if a is in the core!
  if (a_in_coreQ) {
    double x_tjap1 = (twoj_a + 1); // no occ_frac here
    int kmax = twoj_a;
    const auto &vaak = m_cint.get_y_ijk(phi_a, phi_a);
    auto irmax = phi_a.pinf;
    const auto &L_ab_k = m_cint.get_angular_L_kiakib_k(ki_a, ki_a);
    for (int k = 0; k <= kmax; k++) {
      if (L_ab_k[k] == 0)
        continue;
      for (std::size_t i = 0; i < irmax; i++) {
        // nb: need to 'cut' here, or fails w/ f states...
        vex_a[i] += -1 * L_ab_k[k] * vaak[k][i] * x_tjap1;
      }
    } // k
  }   // if a in core
}

//******************************************************************************
const std::vector<double> &HartreeFock::get_vex(const DiracSpinor &psi) const {
  bool valenceQ{};
  auto i = p_wf->getStateIndex(psi.n, psi.k, valenceQ);
  return valenceQ ? appr_vex_val[i] : appr_vex_core[i];
}

//******************************************************************************
DiracSpinor HartreeFock::vex_psia(const DiracSpinor &phi_a) const
// calculates V_ex Psi_a (returns new Dirac Spinor)
// Psi_a can be any orbital (so long as coulomb integrals exist!)
{
  DiracSpinor vexPsi(phi_a.n, phi_a.k, *(phi_a.p_rgrid));
  vex_psia(phi_a, vexPsi);
  return vexPsi;
}
void HartreeFock::vex_psia(const DiracSpinor &phi_a, DiracSpinor &vexPsi) const
// calculates V_ex Psi_a (returns new Dirac Spinor)
// Psi_a can be any orbital (so long as coulomb integrals exist!)
{
  vexPsi *= 0.0;
  if (m_excludeExchange)
    return;

  auto ki_a = phi_a.k_index();
  auto twoj_a = phi_a.twoj();
  std::size_t init = phi_a.k == -1 ? 0 : 1; //?
  for (const auto &phi_b : p_wf->core_orbitals) {
    auto tjb = phi_b.twoj();
    double x_tjbp1 = (tjb + 1) * phi_b.occ_frac;
    auto irmax = std::min(phi_a.pinf, phi_b.pinf);
    int kmin = std::abs(twoj_a - tjb) / 2;
    int kmax = (twoj_a + tjb) / 2;
    const auto &vabk = m_cint.get_y_ijk(phi_b, phi_a);
    const auto &L_ab_k = m_cint.get_angular_L_kiakib_k(ki_a, phi_b.k_index());
    for (int k = kmin; k <= kmax; k++) {
      if (L_ab_k[k - kmin] == 0)
        continue;
      for (auto i = init; i < irmax; i++) {
        auto v = -x_tjbp1 * L_ab_k[k - kmin] * vabk[k - kmin][i];
        vexPsi.f[i] += v * phi_b.f[i];
        vexPsi.g[i] += v * phi_b.g[i];
      } // r
    }   // k
  }     // b
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//******************************************************************************
void HartreeFock::iterate_core_orbital(
    DiracSpinor &phi, const std::vector<double> &vl,
    const DirectHamiltonian &Hd, const std::vector<double> &fVdir0) const {

  const auto eps_target = m_eps_HF;

  const int n = phi.n;
  const int k = phi.k;

  const auto &vd = Hd.getVdir();

  double prev_inner_eps = 1.0;
  auto phi0 = DiracSpinor(n, k, *(phi.p_rgrid));
  auto phiI = DiracSpinor(n, k, *(phi.p_rgrid));
  auto vexPsi = DiracSpinor(n, k, *(phi.p_rgrid));
  for (int ini = 0; ini <= MAX_HART_ITS; ini++) {

    vex_psia(phi, vexPsi);
    const auto Sr = (fVdir0 * phi) - (vd * phi) - vexPsi;

    auto en = Hd.matrixEl(phi, phi) + (phi * vexPsi);
    // if (en > 0) {
    //   // XXX Need for Li? Why??
    //   // XXX Doesn't solve issue..
    //   phi.eps = 1.0;
    //   phi.its = MAX_HART_ITS;
    //   break;
    // }

    auto inner_eps = fabs((en - phi.en) / phi.en);
    bool inner_converged = (inner_eps <= eps_target || ini == MAX_HART_ITS);
    bool inner_getting_worse = (inner_eps > 1.1 * prev_inner_eps && ini > 10);

    if (prev_inner_eps > inner_eps)
      prev_inner_eps = inner_eps;
    if (inner_converged || inner_getting_worse) {
      phi.en = en;
      phi.eps = inner_eps;
      phi.its = ini;
      break;
    }

    solve_inhomog_Green(phi, phi0, phiI, en, vl, p_wf->get_alpha(), Sr);
    phi.normalise();
    phi.en = en;
  }
  DEBUG(printf(" refine--- %2i,%2i: en=%11.5f  HFeps = %.0e;  %2i  "
               "(%4i)\n",
               phi.n, phi.k, phi.en, phi.eps, phi.its, (int)phi.pinf);)
}

//******************************************************************************
inline void HartreeFock::refine_valence_orbital_exchange(DiracSpinor &phi) {

  auto eps_target = m_eps_HF;
  const auto a_damp = 0.1;

  const int n = phi.n;
  const int k = phi.k;

  auto vl = p_wf->vnuc;
  for (unsigned i = 0; i < vl.size(); i++) {
    vl[i] += p_wf->vdir[i];
  }

  DirectHamiltonian Hd(p_wf->vnuc, p_wf->vdir, p_wf->get_alpha());

  auto prev_en = phi.en;
  m_cint.form_core_valence(phi); // only needed if not already done!
  double eps_prev = 1.0;
  auto en = phi.en;
  auto phi0 = DiracSpinor(n, k, p_wf->rgrid);
  auto phiI = DiracSpinor(n, k, p_wf->rgrid);
  auto vexPsi = DiracSpinor(n, k, p_wf->rgrid);
  for (int it = 0; it <= MAX_HART_ITS; ++it) {
    vex_psia(phi, vexPsi);
    auto Sr = -1.0 * vexPsi;
    en = Hd.matrixEl(phi, phi) + phi * vexPsi;

    auto eps = fabs((prev_en - en) / en);
    bool getting_worse = (it > 20 && eps >= 2.0 * eps_prev);
    bool converged = (eps <= eps_target && it > 0);
    if (converged || getting_worse || it == MAX_HART_ITS) {
      phi.en = en;
      phi.eps = eps;
      phi.its = it;
      if (verbose)
        printf("refine: %2i %2i | %3i eps=%6.1e  en=%11.8f\n", n, k, it, eps,
               phi.en);
      break;
    }
    if (eps < eps_prev)
      eps_prev = eps;
    prev_en = en;

    auto oldphi = phi;
    solve_inhomog_Green(phi, phi0, phiI, en, vl, p_wf->get_alpha(), Sr);
    phi.normalise();
    phi = a_damp * oldphi + (1.0 - a_damp) * phi;
    phi.en = en;
    p_wf->orthonormaliseWrtCore(phi);
    m_cint.form_core_valence(phi);
  }
  p_wf->orthonormaliseWrtCore(phi);
  return;
}

//******************************************************************************
inline void HartreeFock::refine_core_orbitals_exchange() {

  const double eps_target = m_eps_HF;

  m_cint.form_core_core(); // only needed if not already done!

  DirectHamiltonian Hd(p_wf->vnuc, p_wf->vdir, p_wf->get_alpha());

  const auto a_damp = 0.2;
  double prev_eps = 1.0;

  for (int it = 0; it <= MAX_HART_ITS; it++) {

    const auto f_core = double(p_wf->Ncore() - 1) / double(p_wf->Ncore());

    auto vl = p_wf->vnuc;
    auto size = vl.size();
    std::vector<double> fVdir0; // allocate once outside?
    fVdir0.reserve(size);
    for (unsigned i = 0; i < size; i++) {
      auto fvd = f_core * p_wf->vdir[i];
      vl[i] += fvd;
      fVdir0.push_back(fvd);
    }
    // const auto vl = vl_tmp;

    std::vector<double> prior_en_list;
    for (auto &phi : p_wf->core_orbitals) {
      // note: can't parallelise -- not thread safe
      prior_en_list.push_back(phi.en);
      auto oldphi = phi; // re-use? nb: n,k each time!
      iterate_core_orbital(phi, vl, Hd, fVdir0);
      phi = (1.0 - a_damp) * phi + a_damp * oldphi;
    }

    double eps = 0.0;
    for (std::size_t i = 0; i < p_wf->core_orbitals.size(); i++) {
      auto en_this = p_wf->core_orbitals[i].en;
      auto en_prev = prior_en_list[i];
      auto d_eps = fabs((en_prev - en_this) / en_this);
      if (d_eps > eps)
        eps = d_eps;
    }
    if (eps < prev_eps)
      prev_eps = eps;

    bool getting_worse = (it > 20 && eps > 1.1 * prev_eps);
    bool converged = ((eps <= eps_target && it > 0) || it == MAX_HART_ITS);
    if (converged || getting_worse) {
      if (verbose)
        printf("refine core  it:%3i eps=%6.1e              \n", it, eps);
      break;
    }
    p_wf->orthonormaliseOrbitals(p_wf->core_orbitals);
    m_cint.form_core_core();
    form_vdir(p_wf->vdir);
    Hd.updateVdir(p_wf->vdir);
  }
  p_wf->orthonormaliseOrbitals(p_wf->core_orbitals, 2);
}

//******************************************************************************
void HartreeFock::solve_inhomog_Green(DiracSpinor &phi, const double en,
                                      const std::vector<double> &v,
                                      const double alpha,
                                      const DiracSpinor &source) const
// make static!?
// NOTE: returns NON-normalised function!
{
  auto phi0 = DiracSpinor(phi.n, phi.k, *phi.p_rgrid);
  auto phiI = DiracSpinor(phi.n, phi.k, *phi.p_rgrid);
  solve_inhomog_Green(phi, phi0, phiI, en, v, alpha, source);
}
//------------------------------------------------------------------------------
void HartreeFock::solve_inhomog_Green(DiracSpinor &phi, DiracSpinor &phi0,
                                      DiracSpinor &phiI, const double en,
                                      const std::vector<double> &v,
                                      const double alpha,
                                      const DiracSpinor &source) const
// make static!?
// NOTE: returns NON-normalised function!
{
  Adams::diracODE_regularAtOrigin(phi0, en, v, alpha);
  Adams::diracODE_regularAtInfinity(phiI, en, v, alpha);
  yfun(phi, phiI, phi0, source);
}

//******************************************************************************
void HartreeFock::yfun(DiracSpinor &phi, const DiracSpinor &phi1,
                       const DiracSpinor &phi2, const DiracSpinor &Sr) const {

  // Wronskian. Note: in current method: only need the SIGN
  int pp = int(0.65 * double(phi1.pinf));
  auto w2 = (phi1.f[pp] * phi2.g[pp] - phi2.f[pp] * phi1.g[pp]);
  // XXX test making this non-constant over r?

  // save typing:
  const auto &gr = *phi.p_rgrid;
  constexpr auto ztr = NumCalc::zero_to_r;
  constexpr auto rti = NumCalc::r_to_inf;

  phi *= 0.0;
  NumCalc::additivePIntegral<ztr>(phi.f, phi1.f, phi2.f, Sr.f, gr, phi1.pinf);
  NumCalc::additivePIntegral<ztr>(phi.f, phi1.f, phi2.g, Sr.g, gr, phi1.pinf);
  NumCalc::additivePIntegral<rti>(phi.f, phi2.f, phi1.f, Sr.f, gr, phi1.pinf);
  NumCalc::additivePIntegral<rti>(phi.f, phi2.f, phi1.g, Sr.g, gr, phi1.pinf);
  NumCalc::additivePIntegral<ztr>(phi.g, phi1.g, phi2.f, Sr.f, gr, phi1.pinf);
  NumCalc::additivePIntegral<ztr>(phi.g, phi1.g, phi2.g, Sr.g, gr, phi1.pinf);
  NumCalc::additivePIntegral<rti>(phi.g, phi2.g, phi1.f, Sr.f, gr, phi1.pinf);
  NumCalc::additivePIntegral<rti>(phi.g, phi2.g, phi1.g, Sr.g, gr, phi1.pinf);
  phi *= (1.0 / w2);
}
