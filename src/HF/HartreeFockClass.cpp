#include "HF/HartreeFockClass.hpp"
#include "Adams/Adams_Greens.hpp"
#include "Adams/DiracODE.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "Dirac/Wavefunction.hpp"
#include "HF/CoulombIntegrals.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/Parametric_potentials.hpp"
#include "Physics/Wigner_369j.hpp"
#include <cmath>
#include <functional>
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
// For non-constant damping
// Slowly ramps the damping factor from a_beg to a_end over interval (beg, end)
static inline auto rampedDamp(double a_beg, double a_end, int beg, int end) {
  return [=](int i) {
    if (i >= end)
      return a_end;
    if (i <= beg)
      return a_beg;
    return (a_end * (i - beg) + a_beg * (end - i)) / (end - beg);
  };
}

//******************************************************************************
HartreeFock::HartreeFock(HFMethod method, Wavefunction &wf,
                         const std::string &in_core, double eps_HF, double h_d,
                         double g_t) //
    : p_wf(&wf),                     //
      p_rgrid(&wf.rgrid),
      m_cint(Coulomb(wf.rgrid, wf.core_orbitals, wf.valence_orbitals)),
      m_eps_HF([=]() { // can give as log..
        return (std::fabs(eps_HF) < 1) ? eps_HF : std::pow(10, -1 * eps_HF);
      }()), //
      m_excludeExchange([=]() {
        return (method == HFMethod::HartreeFock || method == HFMethod::ApproxHF)
                   ? false
                   : true;
      }()),            //
      m_method(method) //
{

  auto param = (method == HFMethod::GreenPRM || method == HFMethod::TietzPRM);
  auto log_eps_init = param ? 15 : 5;

  if (wf.core_orbitals.empty() || param)
    starting_approx_core(in_core, log_eps_init, method, h_d, g_t);

  m_cint.initialise_core_core();
  appr_vex_core.resize(p_wf->core_orbitals.size(),
                       std::vector<double>(p_rgrid->ngp));

  // XXX Update this so that you can create HF class, THEN solve for core later
  switch (method) {
  case HFMethod::HartreeFock:
    hf_core_approx(1.0e-5);
    hf_core_refine();
    break;
  case HFMethod::ApproxHF:
    hf_core_approx(m_eps_HF);
    break;
  case HFMethod::Hartree:
    hf_core_approx(m_eps_HF);
    break;
  default:
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
HFMethod HartreeFock::parseMethod(const std::string &in_method) {
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
void HartreeFock::starting_approx_core(const std::string &in_core,
                                       int log_converge, HFMethod method,
                                       double h_g, double d_t)
// Starting approx for HF. Uses Green parametric
// Later, can put other options if you want.
{
  if (method != HFMethod::TietzPRM) {
    if (std::fabs(h_g * h_g) < 1.0e-6)
      Parametric::defaultGreenCore(p_wf->Znuc(), h_g, d_t);
    p_wf->vdir = Parametric::GreenPotential(p_wf->Znuc(), p_rgrid->r, h_g, d_t);
  } else {
    if (std::fabs(h_g * h_g) < 1.0e-6)
      Parametric::defaultTietz(p_wf->Znuc(), d_t, h_g); // nb: g/t swaped!
    p_wf->vdir = Parametric::TietzPotential(p_wf->Znuc(), p_rgrid->r, h_g, d_t);
  }

  p_wf->solveInitialCore(in_core, log_converge);
  if (p_wf->core_orbitals.empty()) {
    // If H-like, kill "initial" vdir (Green potential)
    p_wf->vdir = std::vector<double>(p_wf->rgrid.ngp, 0);
  }
}

//******************************************************************************
EpsIts HartreeFock::hf_valence(DiracSpinor &phi, std::vector<double> &vexa) {

  auto do_refine =
      (m_method == HFMethod::HartreeFock && !p_wf->core_orbitals.empty());
  auto eps_target_HF = do_refine ? 1.0e-5 : m_eps_HF;

  auto ei = hf_valence_approx(phi, vexa, eps_target_HF);
  if (do_refine)
    ei = hf_valence_refine(phi);
  return ei;
}

//******************************************************************************
void HartreeFock::hf_core_approx(const double eps_target_HF) {

  if (p_wf->core_orbitals.empty()) {
    return;
  }

  auto damper = rampedDamp(0.7, 0.2, 3, 10);
  // don't include all pts in PT for new e guess:
  static const std::size_t de_stride = 5;

  // initialise 'old' potentials
  auto vdir_old = p_wf->vdir;
  auto vex_old = appr_vex_core;

  // Start the HF itterative procedure:
  int hits = 1;
  double t_eps = 1.0;
  auto t_eps_prev = 1.0;
  for (; hits < m_max_hf_its; hits++) {
    DEBUG(std::cerr << "HF core it: " << hits << "\n";)
    auto eta = damper(hits);

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
      p_wf->vdir[j] = (1.0 - eta) * p_wf->vdir[j] + eta * vdir_old[j];
      for (std::size_t i = 0; i < p_wf->core_orbitals.size(); i++) {
        appr_vex_core[i][j] =
            (1.0 - eta) * appr_vex_core[i][j] + eta * vex_old[i][j];
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
      t_eps = (state_eps > t_eps) ? state_eps : t_eps;
      DEBUG(printf(" --- %2i,%2i: en=%11.5f  HFeps = %.0e;  Adams = %.0e[%2i]  "
                   "(%4i)\n",
                   phi.n, phi.k, phi.en, state_eps, phi.eps, phi.its,
                   (int)phi.pinf);)
    } // core states
    DEBUG(std::cerr << "HF core it: " << hits << ": eps=" << t_eps << "\n\n";
          std::cin.get();)

    // Force all core orbitals to be orthogonal to each other
    if (m_explicitOrthog_cc)
      p_wf->orthonormaliseOrbitals(p_wf->core_orbitals, 1);
    auto getting_worse = (hits > 20 && t_eps > t_eps_prev && t_eps < 1.e-5);
    auto converged = (t_eps < eps_target_HF);
    if (converged || getting_worse)
      break;
    t_eps_prev = t_eps;
  } // hits
  if (verbose && m_method != HFMethod::HartreeFock)
    printf("HF core      it:%3i eps=%6.1e\n", hits, t_eps);

  // Now, re-solve core orbitals with higher precission
  for (std::size_t i = 0; i < p_wf->core_orbitals.size(); i++) {
    p_wf->solveDirac(p_wf->core_orbitals[i], p_wf->core_orbitals[i].en,
                     appr_vex_core[i], 15);
  }
  if (m_explicitOrthog_cc)
    p_wf->orthonormaliseOrbitals(p_wf->core_orbitals, 2);
}

//******************************************************************************
void HartreeFock::solveValence() {

  const auto Nval = p_wf->valence_orbitals.size();
  appr_vex_val.resize(Nval);
  m_cint.initialise_core_valence();

  std::vector<EpsIts> eis(Nval);

#pragma omp parallel for
  for (std::size_t i = 0; i < Nval; i++) {
    eis[i] = hf_valence(p_wf->valence_orbitals[i], appr_vex_val[i]);
  }

  double eps_worst = 0.0, eps_best = 10.0;
  std::size_t i_worst = 0, i_best = 0;
  for (std::size_t i = 0; i < Nval; i++) {
    const auto &ei = eis[i];
    if (ei.eps >= eps_worst) {
      eps_worst = ei.eps;
      i_worst = i;
    }
    if (ei.eps < eps_best) {
      eps_best = ei.eps;
      i_best = i;
    }
  }
  if (verbose)
    printf("HF valence:  %3i eps=%6.1e for %s  [%6.1e for %s w/%3i]\n", //
           eis[i_worst].its, eis[i_worst].eps,
           p_wf->valence_orbitals[i_worst].symbol().c_str(), eis[i_best].eps,
           p_wf->valence_orbitals[i_best].symbol().c_str(), eis[i_best].its);
}

//******************************************************************************
EpsIts HartreeFock::hf_valence_approx(DiracSpinor &phi,
                                      std::vector<double> &vexa,
                                      double eps_target_HF)
// Solves HF for given orbital phi, in frozen core.
// Does not store vex (must be done outside)
// Can be used to generate a set of virtual/basis orbitals
{
  phi.occ_frac = 1. / phi.twojp1();

  auto damper = rampedDamp(0.7, 0.2, 2, 6);
  // don't include all pts in PT for new e guess
  static const std::size_t de_stride = 5;

  vexa.clear();
  vexa.resize(p_rgrid->ngp, 0);

  auto vexa_old = vexa;

  p_wf->solveDirac(phi, 0, 15);

  double eps = -1, eps_prev = -1;
  int hits = 1;
  for (; hits < m_max_hf_its; hits++) {
    auto eta = damper(hits);

    double en_old = phi.en;
    vexa_old = vexa;

    m_cint.form_core_valence(phi);
    form_approx_vex_a(phi, vexa);

    for (std::size_t i = 0; i < p_rgrid->ngp; i++) {
      vexa[i] = (1.0 - eta) * vexa[i] + eta * vexa_old[i];
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
    // Force valence state to be orthogonal to core:
    if (m_explicitOrthog_cv)
      p_wf->orthonormaliseWrtCore(phi);

    auto getting_worse = (hits > 20 && eps >= eps_prev && eps < 1.e-5);
    auto converged = (eps <= eps_target_HF);
    if (converged || getting_worse)
      break;
    eps_prev = eps;
  }
  DEBUG(printf("HF val: %2i %2i | %3i eps=%6.1e  en=%11.8f\n", phi.n, phi.k,
               hits, eps, phi.en);)

  if (m_explicitOrthog_cv)
    p_wf->orthonormaliseWrtCore(phi);
  return {eps, hits};
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
  for (const auto &phi_a : p_wf->core_orbitals) {
    const auto tja = phi_a.twoj();

    double e1 = 0, e2 = 0, e3 = 0;
    const double xtjap1 = (tja + 1) * phi_a.occ_frac;
    e1 += xtjap1 * phi_a.en;
    for (const auto &phi_b : p_wf->core_orbitals) {
      const auto tjb = phi_b.twoj();
      const double xtjbp1 = (tjb + 1) * phi_b.occ_frac;
      const auto irmax = std::min(phi_a.pinf, phi_b.pinf);
      const auto &v0bb = m_cint.get_y_ijk(phi_b, phi_b, 0);
      auto R0f2 = NumCalc::integrate(
          {&phi_a.f, &phi_a.f, &v0bb, &p_rgrid->drdu}, 1.0, 0, irmax);
      auto R0g2 = NumCalc::integrate(
          {&phi_a.g, &phi_a.g, &v0bb, &p_rgrid->drdu}, 1.0, 0, irmax);
      e2 += xtjap1 * xtjbp1 * (R0f2 + R0g2);
      // take advantage of symmetry for third term:
      if (phi_b > phi_a)
        continue;
      const double y = (phi_a == phi_b) ? 1.0 : 2.0;
      const int kmin = std::abs(tja - tjb) / 2;
      const int kmax = (tja + tjb) / 2;
      const auto &vabk = m_cint.get_y_ijk(phi_a, phi_b);
      const auto &L_abk =
          m_cint.get_angular_L_kiakib_k(phi_a.k_index(), phi_b.k_index());
      for (auto k = kmin; k <= kmax; k++) {
        if (L_abk[k - kmin] == 0)
          continue;
        const auto ik = k - kmin;
        double R0f3 =
            NumCalc::integrate({&phi_a.f, &phi_b.f, &vabk[ik], &p_rgrid->drdu});
        double R0g3 =
            NumCalc::integrate({&phi_a.g, &phi_b.g, &vabk[ik], &p_rgrid->drdu});
        e3 += y * xtjap1 * xtjbp1 * L_abk[k - kmin] * (R0f3 + R0g3);
      }
    }
    Etot += e1 - 0.5 * (e2 - e3) * p_rgrid->du; // update running total
  }
  return Etot;
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
  const double sf = re_scale ? (1. - 1. / p_wf->Ncore()) : 1;
  for (const auto &phi_b : p_wf->core_orbitals) {
    const double f_sf = sf * (phi_b.twoj() + 1) * phi_b.occ_frac;
    const auto &v0bb = m_cint.get_y_ijk(phi_b, phi_b, 0);
    for (std::size_t i = 0; i < p_rgrid->ngp; i++) {
      vdir[i] += v0bb[i] * f_sf;
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
// psi_a is small. Luckily, however, we don't care as much when psi_a is
// small! Also, since v_ex is already small (compared to vdir), we can make
// good approximation. Therefore, I only calculate vex_a when a=b, or when
// |psi_a| > 1.e3 Further, largest part of v_ex is when a=b. In this case, the
// factor=1 is exact!
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
// calculates V_ex Psi_a
// Psi_a can be any orbital (so long as coulomb integrals exist!)
{
  vexPsi.pinf = phi_a.f.size(); // silly hack. Make sure vexPsi = 0 after pinf
  vexPsi *= 0.0;
  vexPsi.pinf = phi_a.pinf;

  if (m_excludeExchange)
    return;

  auto ki_a = phi_a.k_index();
  auto twoj_a = phi_a.twoj();
  std::size_t init = 1; // phi_a.k == -1 ? 0 : 1; //? Fails f-states otherwise??
  for (const auto &phi_b : p_wf->core_orbitals) {
    auto tjb = phi_b.twoj();
    double x_tjbp1 = (phi_a == phi_b) ? (tjb + 1) : (tjb + 1) * phi_b.occ_frac;
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

// -----------------------------------------------------------------------------
DiracSpinor HartreeFock::vex_psia_any(const DiracSpinor &phi_a,
                                      const std::vector<DiracSpinor> &core,
                                      int k_cut) // static
// calculates V_ex Psi_a (returns new Dirac Spinor)
// Psi_a can be any orbital (Calculates coulomb integrals here!)
{
  DiracSpinor vexPsi(phi_a.n, phi_a.k, *(phi_a.p_rgrid));
  vexPsi.pinf = phi_a.pinf;

  std::vector<double> vabk(phi_a.p_rgrid->ngp);
  // XXX ALSO move this!

  auto tja = phi_a.twoj();
  auto la = phi_a.l();
  std::size_t init = 1; //?
  for (const auto &phi_b : core) {
    auto tjb = phi_b.twoj();
    auto lb = phi_b.l();
    double x_tjbp1 = (phi_a == phi_b) ? (tjb + 1) : (tjb + 1) * phi_b.occ_frac;
    auto irmax = std::min(phi_a.pinf, phi_b.pinf);
    int kmin = std::abs(tja - tjb) / 2;
    if (kmin > k_cut)
      continue;
    int kmax = (tja + tjb) / 2;
    if (kmax > k_cut)
      kmax = k_cut;
    for (int k = kmin; k <= kmax; k++) {
      auto parity = Wigner::parity(la, lb, k);
      if (parity == 0)
        continue;
      auto tjs = Wigner::threej_2(tjb, tja, 2 * k, -1, 1, 0); // XXX lookup!
      if (tjs == 0)
        continue;
      Coulomb::calculate_y_ijk(phi_b, phi_a, k, vabk);
      for (auto i = init; i < irmax; i++) {
        auto v = -x_tjbp1 * tjs * tjs * vabk[i];
        vexPsi.f[i] += v * phi_b.f[i];
        vexPsi.g[i] += v * phi_b.g[i];
      } // r
    }   // k
  }     // b

  return vexPsi;
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//******************************************************************************
void HartreeFock::hf_orbital(DiracSpinor &phi, double en,
                             const std::vector<double> &vl,
                             const DiracSpinor &vx_phi,
                             const std::vector<DiracSpinor> &static_core,
                             const std::vector<double> &v0) const
// Solve Dirac Equation (Eigenvalue): (move to DiracODE??)
//  (H0 + Vl + Vx)Psi = 0
//  (H0 + Vl)Psi = -VxPsi
// Vl is local (e.g., Vnuc + fVdir), Vx is non-local (e.g., (1-f)Vdir + Vex)
// where v0 = (1-f)Vdir  [f=1 for valence states!, so v0 may be empty]
// Small energy adjustmenets (and wfs), solve:
// (Hl - e) dF = de * F -VxPsi
// e -> e+de, F->F+dF
// Core is input so can call in a thread-safe way! (with a 'old_core' copy)
// Only used in dE from dF
{

  // pull these outside? But make sure thread safe!
  DiracSpinor phi0(phi.n, phi.k, *(phi.p_rgrid));
  DiracSpinor phiI(phi.n, phi.k, *(phi.p_rgrid));
  DiracSpinor VxFh(phi.n, phi.k, *(phi.p_rgrid));

  const auto alpha = p_wf->get_alpha();
  DiracODE::solve_inhomog(phi, phi0, phiI, en, vl, alpha, -1.0 * vx_phi);

  // make small adjustments to energy to normalise psi:
  DiracSpinor del_phi(phi.n, phi.k, *(phi.p_rgrid));
  DiracODE::Adams::GreenSolution(del_phi, phiI, phi0, alpha, phi);
  // should del_phi = del_e * del_phi, but makes it worse?
  // nb: after first it, becomes correct.
  auto del_E = 0.5 * (phi * phi - 1.0) / (phi * del_phi);
  auto eps = std::abs(del_E / en);
  int tries = 0;
  for (; tries < m_max_hf_its; ++tries) { // m_max_hf_its
    if (eps < m_eps_HF)
      break;
    {
      if (!v0.empty()) {     // essentially, for core:
        VxFh = v0 * del_phi; // v0 = (1-f)Vd;
      } else {               // essentially, for valence
        VxFh = vex_psia_any(del_phi, static_core, 0);
      }
      // Depends very strengely on VxFh...
      DiracODE::Adams::GreenSolution(del_phi, phiI, phi0, alpha,
                                     del_E * phi - VxFh);
    }
    auto delta_Norm = phi * phi - 1.0;
    auto de0 = del_E;
    del_E = 0.5 * de0 * delta_Norm / (phi * del_phi);
    eps = std::abs(del_E / en);
    en += del_E;
    phi -= (1.0 * del_E / de0) * del_phi;
  }
  phi.en = en;
  phi.eps = eps;
  phi.its = tries;
  // phi.normalise(); //? Not needed
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//******************************************************************************
EpsIts HartreeFock::hf_valence_refine(DiracSpinor &phi) {

  if (p_wf->core_orbitals.empty())
    return {0, 0};

  const auto eps_target = m_eps_HF;

  auto damper = rampedDamp(0.7, 0.1, 5, 25);
  double extra_damp = 0.0;

  const auto vl = NumCalc::sumVecs({&(p_wf->vnuc), &(p_wf->vdir)});

  const auto phi_zero = phi;
  const auto vexPsi_zero = get_vex(phi) * phi;

  auto prev_en = phi.en;
  m_cint.form_core_valence(phi); // only needed if not already done!
  double best_eps = 1.0;
  auto phi0 = DiracSpinor(phi.n, phi.k, p_wf->rgrid);
  auto phiI = DiracSpinor(phi.n, phi.k, p_wf->rgrid);
  auto vexPsi = DiracSpinor(phi.n, phi.k, p_wf->rgrid);
  int it = 0;
  double eps = 1.0;
  int worse_count = 0;
  for (; it <= m_max_hf_its; ++it) {
    auto a_damp = damper(it) + extra_damp;

    vex_psia(phi, vexPsi);
    auto oldphi = phi;
    auto en = phi_zero.en +
              (phi_zero * vexPsi - phi * vexPsi_zero) / (phi * phi_zero);
    hf_orbital(phi, en, vl, vexPsi, p_wf->core_orbitals);
    eps = std::fabs((prev_en - phi.en) / phi.en);
    prev_en = phi.en;

    if (it > 20 && eps > 1.5 * best_eps) {
      ++worse_count;
      extra_damp = extra_damp > 0 ? 0 : 0.1;
    } else {
      worse_count = 0;
    }
    const bool converged = (eps <= eps_target && it > 0);
    if (converged || worse_count > 2)
      break;

    if (eps < best_eps)
      best_eps = eps;

    DEBUG(std::cout << it << " " << eps << " " << phi.en << " "
                    << en - phi_zero.en << " " << phi * phi << "\n";)

    phi = (1.0 - a_damp) * phi + a_damp * oldphi;
    if (m_explicitOrthog_cv) {
      p_wf->orthonormaliseWrtCore(phi);
    } else {
      phi.normalise();
    }
    m_cint.form_core_valence(phi);
  } // End HF its

  if (m_explicitOrthog_cv)
    p_wf->orthonormaliseWrtCore(phi);

  DEBUG(printf("refine: %2i %2i | %3i eps=%6.1e  en=%11.8f\n", phi.n, phi.k, it,
               eps, phi.en);)
  return {eps, it};
}

//******************************************************************************
inline void HartreeFock::hf_core_refine() {

  if (p_wf->core_orbitals.empty()) {
    return;
  }

  const double eps_target = m_eps_HF;
  m_cint.form_core_core(); // only needed if not already done!
  auto damper = rampedDamp(0.8, 0.1, 5, 30);
  double extra_damp = 0;

  std::vector<double> vl(p_wf->rgrid.ngp); // Vnuc + fVd
  std::vector<double> v0(p_wf->rgrid.ngp); // (1-f)Vd
  const auto f_core = double(p_wf->Ncore() - 1) / double(p_wf->Ncore());
  const auto &vd = p_wf->vdir;

  // Store arrays of intitial Psi and VexPsi, and VdirPsi (for En guess)
  // And allocate arrays for VexPsi, so can //-ise it loop (over orbs)!
  const auto core_zero = p_wf->core_orbitals;
  auto core_prev = p_wf->core_orbitals;
  std::vector<DiracSpinor> vexCore_zero;
  const auto vd0 = p_wf->vdir;
  std::vector<DiracSpinor> vexF_list;
  const auto Ncore = p_wf->core_orbitals.size();
  std::vector<double> eps_lst(Ncore, 0.0);
  for (std::size_t i = 0; i < Ncore; ++i) {
    auto &phi = p_wf->core_orbitals[i];
    vexCore_zero.push_back(get_vex(phi) * phi);
    vexF_list.push_back(DiracSpinor(phi.n, phi.k, *(phi.p_rgrid)));
  }

  double eps = 0.0;
  double best_eps = 1.0;
  double best_worst_eps = 1.0;
  std::size_t worst_index = 0;
  std::size_t best_index = 0;
  int worse_count = 0;
  int it = 0;
  for (; it <= m_max_hf_its; it++) {
    auto a_damp = damper(it) + extra_damp;

    // re-calculate each Vl = vnuc + fvdir, v0 = (1-f)vdir:
    for (auto i = 0ul; i < p_wf->rgrid.ngp; i++) {
      vl[i] = p_wf->vnuc[i] + f_core * vd[i];
      v0[i] = (1.0 - f_core) * vd[i];
    }

    // re-calculate each VexPsi:
    for (std::size_t i = 0; i < Ncore; ++i) {
      const auto &phi = p_wf->core_orbitals[i];
      vex_psia(phi, vexF_list[i]);
    }

    core_prev = p_wf->core_orbitals;

#pragma omp parallel for
    for (std::size_t i = 0; i < Ncore; ++i) {
      auto &phi = p_wf->core_orbitals[i];
      const auto &phi_zero = core_zero[i];
      const auto &vexPsi_zero = vexCore_zero[i];

      // const auto oldphi = phi;
      const auto oldphi = core_prev[i];
      const auto &vexPsi = vexF_list[i];
      auto en = phi_zero.en + (phi_zero * vexPsi - phi * vexPsi_zero +
                               phi_zero * (vd * phi) - phi * (vd0 * phi_zero)) /
                                  (phi * phi_zero);
      const auto v_nonlocal = v0 * phi + vexPsi;
      hf_orbital(phi, en, vl, v_nonlocal, core_prev, v0);
      phi = (1.0 - a_damp) * phi + a_damp * oldphi;
      phi.normalise();
      auto d_eps = std::fabs((oldphi.en - phi.en) / phi.en);
      eps_lst[i] = d_eps;
    }

    eps = eps_lst[0];
    best_eps = eps_lst[0];
    for (std::size_t i = 1; i < Ncore; ++i) {
      auto t_eps = eps_lst[i];
      if (t_eps >= eps) {
        eps = t_eps;
        worst_index = i;
      }
      if (t_eps < best_eps) {
        best_eps = t_eps;
        best_index = i;
      }
    }
    DEBUG(std::cout << eps_target << " " << eps << " "
                    << p_wf->core_orbitals[worst_index].symbol() << " -- "
                    << " " << beps << " "
                    << p_wf->core_orbitals[bindex].symbol() << "\n";)

    if (it > 20 && eps > 1.5 * best_worst_eps) {
      ++worse_count;
      extra_damp = extra_damp > 0 ? 0 : 0.4;
    } else {
      worse_count = 0;
    }
    const bool converged = (eps <= eps_target && it > 0);
    if (converged || worse_count > 3)
      break;

    if (eps < best_worst_eps)
      best_worst_eps = eps;
    if (m_explicitOrthog_cc)
      p_wf->orthonormaliseOrbitals(p_wf->core_orbitals);
    m_cint.form_core_core();
    form_vdir(p_wf->vdir);
  }
  if (m_explicitOrthog_cc)
    p_wf->orthonormaliseOrbitals(p_wf->core_orbitals, 2);

  if (verbose)
    printf("HF core:  it:%3i eps=%6.1e for %s  [%6.1e for %s]\n", //
           it, eps, p_wf->core_orbitals[worst_index].symbol().c_str(), best_eps,
           p_wf->core_orbitals[best_index].symbol().c_str());
}
