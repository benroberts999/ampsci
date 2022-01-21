#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/HartreeFock.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

//******************************************************************************
ContinuumOrbitals::ContinuumOrbitals(const Wavefunction &wf, int izion)
    : rgrid(wf.rgrid),
      p_hf(wf.getHF()),
      Z(wf.Znuc()),
      Zion(izion),
      alpha(wf.alpha),
      v_local(qip::add(wf.vnuc, wf.vdir)) {}

//******************************************************************************
double ContinuumOrbitals::check_orthog(bool print) const {
  double worst = 0.0;
  if (p_hf == nullptr)
    return worst;
  for (const auto &Fc : orbitals) {
    for (const auto &Fn : p_hf->get_core()) {
      if (Fn.k != Fc.k)
        continue;
      const auto eps = Fc * Fn;
      if (std::abs(eps) > std::abs(worst))
        worst = eps;
      if (print) {
        std::cout << "<" << Fc.shortSymbol() << "|" << Fn.shortSymbol()
                  << "> = ";
        printf("%.1e\n", eps);
      }
    }
  }
  return worst;
}

//******************************************************************************
int ContinuumOrbitals::solveContinuumHF(double ec, int max_l,
                                        bool force_rescale, bool subtract_self,
                                        bool force_orthog,
                                        const DiracSpinor *p_psi)
// Overloaded, assumes min_l=0
{
  return solveContinuumHF(ec, 0, max_l, force_rescale, subtract_self,
                          force_orthog, p_psi);
}

//******************************************************************************
int ContinuumOrbitals::solveContinuumHF(double ec, int min_l, int max_l,
                                        bool force_rescale, bool subtract_self,
                                        bool force_orthog,
                                        const DiracSpinor *p_psi)
// Solved the Dirac equation for local potential for positive energy (no mc2)
// continuum (un-bound) states [partial waves].
//  * Goes well past num_points, looks for asymptotic region, where wf is
//  sinosoidal
//  * Uses fit to known exact H-like for normalisation.
{

  // Find 'inital guess' for asymptotic region:
  const double lam = 1.0e6;
  const double r_asym =
      (Zion + std::sqrt(4.0 * lam * ec + std::pow(Zion, 2))) / (2.0 * ec);

  // Check if 'h' is small enough for oscillating region:
  const double h_target = (M_PI / 15) / std::sqrt(2.0 * ec);
  const auto h = rgrid->du();
  if (h > h_target) {
    std::cout << "WARNING 61 CntOrb: Grid not dense enough for ec=" << ec
              << " (du=" << h << ", need du<" << h_target << ")\n";
    if (h > 2 * h_target) {
      std::cout << "FAILURE 64 CntOrb: Grid not dense enough for ec=" << ec
                << " (du=" << h << ", need du<" << h_target << ")\n";
      return 1;
    }
  }

  // nb: Don't need to extend grid each time... but want thread-safe
  auto cgrid = *rgrid;
  cgrid.extend_to(1.1 * r_asym);

  auto vc = v_local;

  // Highest core state (later: use the state being ionised)
  // const auto &Fa = p_hf->get_core().back();

  // Subtracting single electron contribution from direct potential
  if ((p_psi != nullptr) && (subtract_self)) {
    // if (subtract_self) {
    std::vector<double> vd_single = Coulomb::yk_ab(*p_psi, *p_psi, 0);
    // v_local = qip::add(v_local, qip::scale(vd_single, -1.0));
    vc = qip::compose(std::minus{}, vc, vd_single);
    std::cout << "Subtracting self-interaction..." << std::endl;
  }

  // "Z_ion" - "actual" (excluding exchange.....)
  const auto z_tmp = std::abs(vc.back() * rgrid->r().back());

  // Extend local (Vnuc+Vdir) potential to new grid
  vc.reserve(cgrid.num_points());
  for (auto i = rgrid->num_points(); i < cgrid.num_points(); i++) {
    vc.push_back(-z_tmp / cgrid.r(i));
  }

  // Re-scale large-r part of local potential, so goes like -1/r large r
  // Note: doesn't inclue exchange..
  // This also kills orthogonality for HF...
  if (force_rescale) {
    // nb: this agrees best with Dzuba, but bad for orthog
    for (auto i = cgrid.num_points() - 1; i != 0; i--) {
      if (vc[i] > -Zion / cgrid.r(i)) {
        vc[i] = -Zion / cgrid.r(i);
      }
    }
  }

  // Technically, eveything above this needs to happen only once...
  // However, the below code takes ~10x longer than this, so doesn't matter much
  //*******************************

  // loop through each kappa state
  for (int k_i = 0; true; ++k_i) {
    const auto kappa = AtomData::kappaFromIndex(k_i);
    const auto l = AtomData::l_k(kappa);
    if (l < min_l)
      continue;
    if (l > max_l)
      break;

    auto &Fc = orbitals.emplace_back(0, kappa, rgrid);
    Fc.set_en() = ec;
    // solve initial, without exchange term
    DiracODE::solveContinuum(Fc, ec, vc, cgrid, r_asym, alpha);

    // Include exchange (Hartree Fock)
    const int max_its = 20;
    const double conv_target = 1.0e-4;
    if (p_hf != nullptr && !p_hf->excludeExchangeQ()) {
      for (int it = 0; it <= max_its; ++it) {
        const auto vx0 = HF::vex_approx(Fc, p_hf->get_core());
        const auto vl = qip::add(vc, vx0);

        // Copy old solution (needed by DiracODE)
        const auto Fc0 = Fc;
        if (p_hf->method() == HF::Method::HartreeFock) {
          auto VxFc = HF::vexFa(Fc, p_hf->get_core()) - vx0 * Fc;
          // Extend onto larger grid
          VxFc.set_f().resize(vc.size());
          VxFc.set_g().resize(vc.size());
          DiracODE::solveContinuum(Fc, ec, vl, cgrid, r_asym, alpha, &VxFc,
                                   &Fc0);
        } else { // HF::Method::ApproxHF)
          DiracODE::solveContinuum(Fc, ec, vl, cgrid, r_asym, alpha);
        }
        // check convergance:
        const auto eps = ((Fc0 - Fc) * (Fc0 - Fc)) / (Fc * Fc);
        if (eps < conv_target || it == max_its) {
          break;
        }
        // Damp the orbital
        Fc = 0.5 * (Fc + Fc0);
      } // it
    }   // if HF

    // Forcing orthogonality between continuum and core states
    if (force_orthog) {
      for (const auto &Fn : p_hf->get_core()) {
        if (Fc.k == Fn.k) {
          std::cout << "Before forcing orthog: <" << Fc.shortSymbol() << "|"
                    << Fn.shortSymbol() << "> = ";
          printf("%.1e\n", Fc * Fn);
          Fc -= (Fn * Fc) * Fn;
          std::cout << "After: <" << Fc.shortSymbol() << "|" << Fn.shortSymbol()
                    << "> = ";
          printf("%.1e\n", Fc * Fn);
          std::cout << std::endl;
        } // if
      }   // loop through core
    }     // if
    // func already exists for this in Wavefunction.cpp
    // const auto &Fn = p_hf->get_core();
    // Wavefunction::orthogonaliseWrt(Fc, Fn);

  } // kappa

  return 0;
}

//******************************************************************************
int ContinuumOrbitals::solveContinuumZeff(double ec, int min_l, int max_l,
                                          double e_core, double n_core,
                                          bool force_orthog)
// Solves Dirac equation for H-like potential (Zeff model)
// Same Zeff as used by DarkARC (eqn B35 of arxiv:1912.08204):
// Zeff = sqrt{I_{njl} eV / 13.6 eV} * n
// au: Zeff = sqrt{2 * I_{njl}} * n
{

  // Find 'inital guess' for asymptotic region:
  const double lam = 1.0e6;
  const double r_asym =
      (Zion + std::sqrt(4.0 * lam * ec + std::pow(Zion, 2))) / (2.0 * ec);

  // Check if 'h' is small enough for oscillating region:
  const double h_target = (M_PI / 15) / std::sqrt(2.0 * ec);
  const auto h = rgrid->du();
  if (h > h_target) {
    std::cout << "WARNING 61 CntOrb: Grid not dense enough for ec=" << ec
              << " (du=" << h << ", need du<" << h_target << ")\n";
    if (h > 2 * h_target) {
      std::cout << "FAILURE 64 CntOrb: Grid not dense enough for ec=" << ec
                << " (du=" << h << ", need du<" << h_target << ")\n";
      return 1;
    }
  }

  // nb: Don't need to extend grid each time... but want thread-safe
  auto cgrid = *rgrid;
  cgrid.extend_to(1.1 * r_asym);

  const double Zeff = std::sqrt(-2.0 * e_core) * n_core;
  std::cout << "Zeff = " << Zeff << " (n = " << n_core << ", e = " << e_core
            << ")" << std::endl;

  std::vector<double> vc;
  vc.resize(cgrid.num_points());
  for (auto i = 0ul; i < cgrid.num_points(); i++) {
    vc[i] = -Zeff / cgrid.r(i);
  }

  // loop through each kappa state
  for (int k_i = 0; true; ++k_i) {
    const auto kappa = AtomData::kappaFromIndex(k_i);
    const auto l = AtomData::l_k(kappa);
    if (l < min_l)
      continue;
    if (l > max_l)
      break;

    auto &Fc = orbitals.emplace_back(0, kappa, rgrid);
    Fc.set_en() = ec;
    // solve initial, without exchange term
    DiracODE::solveContinuum(Fc, ec, vc, cgrid, r_asym, alpha);

    // Forcing orthogonality between continuum and core states
    if (force_orthog) {
      for (const auto &Fn : p_hf->get_core()) {
        if (Fc.k == Fn.k) {
          std::cout << "Before forcing orthog: <" << Fc.shortSymbol() << "|"
                    << Fn.shortSymbol() << "> = ";
          printf("%.1e\n", Fc * Fn);
          Fc -= (Fn * Fc) * Fn;
          std::cout << "After: <" << Fc.shortSymbol() << "|" << Fn.shortSymbol()
                    << "> = ";
          printf("%.1e\n", Fc * Fn);
          std::cout << std::endl;
        } // if
      }   // loop through core
    }     // if
  }

  return 0;
}

//******************************************************************************
void ContinuumOrbitals::clear() { orbitals.clear(); }
