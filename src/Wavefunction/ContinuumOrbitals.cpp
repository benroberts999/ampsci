#include "Wavefunction/ContinuumOrbitals.hpp"
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
int ContinuumOrbitals::solveContinuumHF(double ec, int max_l)
// Overloaded, assumes min_l=0
{
  return solveContinuumHF(ec, 0, max_l);
}

//******************************************************************************
int ContinuumOrbitals::solveContinuumHF(double ec, int min_l, int max_l)
// Solved the Dirac equation for local potential for positive energy (no mc2)
// continuum (un-bound) states [partial waves].
//  * Goes well past num_points, looks for asymptotic region, where wf is
//  sinosoidal
//  * Uses fit to known exact H-like for normalisation.
{

  // Find 'inital guess' for asymptotic region:
  const double lam = 1.0e7;
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

  // XXX Don't need to extend grid each time...
  // ExtendedGrid cgrid(*rgrid, 1.2 * r_asym);
  auto cgrid = *rgrid;
  cgrid.extend_to(1.2 * r_asym);

  // "Z_ion" - "actual" (excluding exchange.....)
  auto z_tmp = std::abs(v_local.back() * rgrid->r().back());
  std::cout << "z_tmp=" << z_tmp << "\n";
  // If ztm is 0, means neutral atom. Effective charge should be 1
  // Exchange doesn't go further than core...
  // This doesn't seem to have any impact, so unimportant
  if (z_tmp < 1)
    z_tmp = 1;

  // Extend local (Vnuc+Vdir) potential to new grid
  auto vc = v_local;
  vc.reserve(cgrid.num_points());
  for (auto i = rgrid->num_points(); i < cgrid.num_points(); i++) {
    vc.push_back(-z_tmp / cgrid.r(i));
  }

  // Re-scale large-r part of local potential, so goes like -1/r large r
  // Note: doesn't inclue exchange..
  // This also kills orthogonality for HF...
  if (force_rescale) {
    // nb: this, without the 'break' agrees best with Dzuba, but bad for orthog
    for (auto i = cgrid.num_points() - 1; i != 0; i--) {
      if (vc[i] > -Zion / cgrid.r(i)) {
        vc[i] = -Zion / cgrid.r(i);
      } else {
        // break; ?
      }
    }
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
    DiracODE::solveContinuum(Fc, ec, vc, cgrid, r_asym, alpha);

    // Include exchange (Hartree Fock)
    if (p_hf != nullptr && !p_hf->excludeExchangeQ()) {
      for (int it = 0; it < 100; ++it) {
        const auto vx0 = HF::vex_approx(Fc, p_hf->get_core());
        const auto vl = qip::add(vc, vx0);
        auto VxFc = HF::vexFa(Fc, p_hf->get_core()) - vx0 * Fc;
        // Extend onto larger grid
        VxFc.set_f().resize(vc.size());
        VxFc.set_g().resize(vc.size());
        // Copy old solution (needed by DiracODE)
        const auto Fc0 = Fc;
        DiracODE::solveContinuum(Fc, ec, vl, cgrid, r_asym, alpha, &VxFc, &Fc0);
        const auto eps = ((Fc0 - Fc) * (Fc0 - Fc)) / (Fc * Fc);
        if (eps < 1.0e-16 || it == 249) {
          // std::cout << Fc.shortSymbol() << " " << it << " " << eps << "\n";
          break;
        }
        Fc = 0.5 * (Fc + Fc0);
      }
    }
  }

  return 0;
}

//******************************************************************************
void ContinuumOrbitals::clear() { orbitals.clear(); }
