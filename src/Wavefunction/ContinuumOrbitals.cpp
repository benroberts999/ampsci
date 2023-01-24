#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/HartreeFock.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

//==============================================================================
ContinuumOrbitals::ContinuumOrbitals(const Wavefunction &wf, int izion)
    : rgrid(wf.grid_sptr()), p_hf(wf.vHF()), Zion(izion), alpha(wf.alpha()) {}

ContinuumOrbitals::ContinuumOrbitals(const HF::HartreeFock *hf, int izion)
    : rgrid(hf->grid_sptr()), p_hf(hf), Zion(izion), alpha(hf->alpha()) {}

//==============================================================================
double ContinuumOrbitals::check_orthog(bool print) const {
  double worst = 0.0;
  if (p_hf == nullptr)
    return worst;
  for (const auto &Fc : orbitals) {
    for (const auto &Fn : p_hf->core()) {
      if (Fn.kappa() != Fc.kappa())
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

//==============================================================================
int ContinuumOrbitals::solveContinuumHF(double ec, int min_l, int max_l,
                                        const DiracSpinor *Fi,
                                        bool force_rescale, bool subtract_self,
                                        bool force_orthog_Fi) {

  // include Hartree here? Probably shouldn't, since we do "core Hartree"
  const auto self_consistant = (p_hf->method() == HF::Method::HartreeFock ||
                                p_hf->method() == HF::Method::ApproxHF);

  // Also orthogonalise against entire core: (makes no difference)
  const bool orthog_core = force_orthog_Fi;

  // Find 'inital guess' for asymptotic region:
  const double r_asym = 2000.0 / std::sqrt(2 * ec);

  // Check if 'du' (step-size at large r) is small enough for oscillating region:
  // nb: Assumes log-linear grid, but should definitely use log-linear grid here!
  const double du_target = (M_PI / 15.0) / std::sqrt(2.0 * ec);
  const auto du = rgrid->du();
  if (du > du_target) {
    std::cout << "WARNING 61 CntOrb: Grid not dense enough for ec=" << ec
              << " (du=" << du << ", need du<" << du_target << ")\n";
    if (du > 2.0 * du_target) {
      std::cout << "FAILURE 64 CntOrb: Grid not dense enough for ec=" << ec
                << " (du=" << du << ", need du<" << du_target << ")\n";
      return 1;
    }
  }

  // nb: Don't need to extend grid each time... but want thread-safe
  auto cgrid = *rgrid;
  cgrid.extend_to(1.1 * r_asym);

  using namespace qip::overloads;
  auto vc = p_hf->vlocal();
  if ((Fi != nullptr) && subtract_self && self_consistant) {
    // Subtract off the self-interaction direct part: V-self(r) = y^0_ii(r)
    vc -= Coulomb::yk_ab(*Fi, *Fi, 0);
  }

  // Find "actual" Z_ion: -V(r)*r ~ Zion/r at large r
  const auto Z_eff = std::max(double(Zion), -vc.back() * rgrid->r().back());

  // Extend local (Vnuc+Vdir-Vself) potential to new grid
  vc.reserve(cgrid.num_points());
  for (auto i = rgrid->num_points(); i < cgrid.num_points(); i++) {
    vc.push_back(-Z_eff / cgrid.r(i));
  }

  // We may wish to do this to test things, but not for final calculations:
  if (force_rescale && subtract_self) {
    std::cout << "\nWarning: should not subtract self interaction _and_ "
                 "rescale V(r): do one or the other\n";
  }

  // Re-scale large-r part of local potential, so goes like -1/r large r
  // Note: doesn't inclue exchange..
  // This also kills orthogonality for HF...
  if (force_rescale && self_consistant) {
    // nb: this agrees best with Dzuba, but bad for orthog
    for (std::size_t i = 0; i < cgrid.num_points(); ++i) {
      if (vc[i] > -Zion / cgrid.r(i)) {
        vc[i] = -Zion / cgrid.r(i);
      }
    }
  }

  // Technically, eveything above this needs to happen only once...
  // However, the below code takes ~10x longer than this, so doesn't matter much
  //==============================*

  // loop through each kappa state
  for (int k_i = 0; true; ++k_i) {
    const auto kappa = Angular::kappaFromIndex(k_i);
    const auto l = Angular::l_k(kappa);
    if (l < min_l)
      continue;
    if (l > max_l)
      break;

    auto &Fc = orbitals.emplace_back(0, kappa, rgrid);
    Fc.en() = ec;
    // solve initial, without exchange term
    DiracODE::solveContinuum(Fc, ec, vc, cgrid, r_asym, alpha);
    // Then, include exchange correction:
    if (p_hf != nullptr && !p_hf->excludeExchangeQ()) {
      IncludeExchange(Fc, Fi, force_orthog_Fi, cgrid, vc, r_asym);
    }

  } // kappa

  // Orthogonalise against entire core:
  if (orthog_core) {
    for (auto &Fc : orbitals) {
      for (const auto &Fa : p_hf->core()) {
        if (Fa.kappa() == Fc.kappa())
          Fc -= (Fc * Fa) * Fa;
      }
      // orthod wrt rest of core can slightly ruin "main" orthog condition
      if (force_orthog_Fi && Fi != nullptr && Fi->kappa() == Fc.kappa()) {
        Fc -= (*Fi * Fc) * *Fi;
      }
    }
  }

  return 0;
}

//******************************************************************************
void ContinuumOrbitals::IncludeExchange(DiracSpinor &Fc, const DiracSpinor *Fi,
                                        bool force_orthog_Fi, const Grid &cgrid,
                                        const std::vector<double> &vc,
                                        double r_asym) {
  // Include exchange (Hartree Fock)
  const int max_its = 50;
  const double conv_target = 1.0e-6;

  for (int it = 0; it <= max_its; ++it) {
    const auto vx0 = HF::vex_approx(Fc, p_hf->core());
    const auto vl = qip::add(vc, vx0);

    // Copy old solution (needed by DiracODE)
    const auto Fc0 = Fc;
    if (p_hf->method() == HF::Method::HartreeFock) {
      auto VxFc = HF::vexFa(Fc, p_hf->core()) - vx0 * Fc;
      // Extend onto larger grid
      VxFc.f().resize(vc.size());
      VxFc.g().resize(vc.size());
      DiracODE::solveContinuum(Fc, Fc.en(), vl, cgrid, r_asym, alpha, &VxFc,
                               &Fc0);
    } else { // HF::Method::ApproxHF)
      DiracODE::solveContinuum(Fc, Fc.en(), vl, cgrid, r_asym, alpha);
    }
    // Force orthogonality to Fi (ionised state) (at each HF step)
    if (force_orthog_Fi && Fi && Fi->kappa() == Fc.kappa()) {
      Fc -= (*Fi * Fc) * *Fi;
    }
    // check convergance:
    const auto eps = ((Fc0 - Fc) * (Fc0 - Fc)) / (Fc * Fc);
    if (eps < conv_target || it == max_its) {
      break;
    }
    // Damp the orbital
    Fc = 0.5 * (Fc + Fc0);
  } // it
}

//******************************************************************************
int ContinuumOrbitals::solveContinuumZeff(double ec, int min_l, int max_l,
                                          double Z_eff, const DiracSpinor *Fi,
                                          bool force_orthog)
// Solves Dirac equation for H-like potential (Zeff model)
// Same Zeff as used by DarkARC (eqn B35 of arxiv:1912.08204):
// Zeff = sqrt{I_{njl} eV / 13.6 eV} * n
// au: Zeff = sqrt{2 * I_{njl}} * n
{

  // Also orthogonalise against entire core: (make no difference)
  const bool orthog_core = force_orthog;

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

  std::vector<double> vc;
  vc.resize(cgrid.num_points());
  for (auto i = 0ul; i < cgrid.num_points(); i++) {
    vc[i] = -Z_eff / cgrid.r(i);
  }

  // loop through each kappa state
  for (int k_i = 0; true; ++k_i) {
    const auto kappa = Angular::kappaFromIndex(k_i);
    const auto l = Angular::l_k(kappa);
    if (l < min_l)
      continue;
    if (l > max_l)
      break;

    auto &Fc = orbitals.emplace_back(0, kappa, rgrid);
    Fc.en() = ec;
    // solve initial, without exchange term
    DiracODE::solveContinuum(Fc, ec, vc, cgrid, r_asym, alpha);

  } // kappa

  // Orthogonalise against entire core:
  if (orthog_core) {
    for (auto &Fc : orbitals) {
      for (const auto &Fa : p_hf->core()) {
        if (Fa.kappa() == Fc.kappa())
          Fc -= (Fc * Fa) * Fa;
      }
    }
  }

  // Forcing orthogonality between continuum states and current core state
  // (have to do this _after_ orthog_core, since that slightly breaks this)
  if (force_orthog) {
    for (auto &Fc : orbitals) {
      if (Fi != nullptr && Fi->kappa() == Fc.kappa()) {
        Fc -= (*Fi * Fc) * *Fi;
      }
    }
  }

  return 0;
}

//==============================================================================
void ContinuumOrbitals::clear() { orbitals.clear(); }
