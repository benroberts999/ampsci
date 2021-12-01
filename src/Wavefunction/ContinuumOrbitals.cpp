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
                                        const DiracSpinor *Fi)
// Overloaded, assumes min_l=0
{
  return solveContinuumHF(ec, 0, max_l, Fi);
}

//******************************************************************************
int ContinuumOrbitals::solveContinuumHF(double ec, int min_l, int max_l,
                                        const DiracSpinor *Fi)
// Solved the Dirac equation for local potential for positive energy (no mc2)
// continuum (un-bound) states [partial waves].
//  * Goes well past num_points, looks for asymptotic region, where wf is
//  sinosoidal
//  * Uses fit to known exact H-like for normalisation.
{

  const bool orthog_Fi = true;
  const bool orthog_core = false;
  const bool subtract_self_int = true;
  const bool force_rescale = false;

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

  // include Hartree here? Probably shouldn't, since we do "core Hartree"
  const auto self_consistant = (p_hf->method() == HF::Method::HartreeFock ||
                                p_hf->method() == HF::Method::ApproxHF
                                /*|| p_hf->method() == HF::Method::Hartree*/
  );

  auto vc = v_local;
  if (Fi && subtract_self_int && self_consistant) {
    // Subtract off the self-interaction direct part
    const auto vdir_sub = Coulomb::yk_ab(*Fi, *Fi, 0);
    qip::compose(std::minus{}, &vc, vdir_sub);
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
        // Orthog (at each HF step)
        if (orthog_Fi && Fi && Fi->k == Fc.k) {
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
    }   // if HF
  }     // kappa

  // Orthogonalise against entire core?
  if (orthog_core) {
    for (auto &phic : orbitals) {
      for (const auto &phi : p_hf->get_core()) {
        if (phic.k == phi.k)
          phic -= (phic * phi) * phi;
      }
    }
  }

  return 0;
}

//******************************************************************************
void ContinuumOrbitals::clear() { orbitals.clear(); }
