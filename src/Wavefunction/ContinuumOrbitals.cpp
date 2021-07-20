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
int ContinuumOrbitals::solveLocalContinuum(double ec, int max_l)
// Overloaded, assumes min_l=0
{
  return solveLocalContinuum(ec, 0, max_l);
}

//******************************************************************************
int ContinuumOrbitals::solveLocalContinuum(double ec, int min_l, int max_l)
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
  const auto z_tmp = std::abs(v_local.back() * rgrid->r().back());

  auto vc = v_local;
  vc.reserve(cgrid.num_points());
  for (auto i = rgrid->num_points(); i < cgrid.num_points(); i++) {
    if (force_rescale) {
      vc.push_back(-Zion / cgrid.r(i));
    } else {
      vc.push_back(-z_tmp / cgrid.r(i));
    }
  }

  for (int k_i = 0; true; ++k_i) { // loop through each kappa state
    const auto kappa = AtomData::kappaFromIndex(k_i);
    const auto l = AtomData::l_k(kappa);
    if (l < min_l)
      continue;
    if (l > max_l)
      break;

    auto &phi = orbitals.emplace_back(0, kappa, rgrid);
    phi.set_en() = ec;
    DiracODE::solveContinuum(phi, ec, vc, cgrid, r_asym, alpha);

    // Include (approximate) exchange
    std::vector<double> vx(vc.size());
    if (!p_hf->excludeExchangeQ()) {
      auto n0 = phi * phi;
      for (int iteration = 0; iteration < 50; ++iteration) {
        qip::add(&vx, HF::vex_approx(phi, p_hf->get_core(), 99, 0.01));
        qip::scale(&vx, 0.5);

        auto vtot = qip::add(vc, vx);

        // Ensure potential goes as - Zion / r at large r
        if (force_rescale) {
          for (std::size_t ir = rgrid->num_points() - 1; ir != 0; --ir) {
            const auto r = rgrid->r()[ir];
            if (r * std::abs(vtot[ir]) > Zion)
              break;
            vtot[ir] = -Zion / r;
          }
        }

        DiracODE::solveContinuum(phi, ec, vtot, cgrid, r_asym, alpha);
        const auto nn = phi * phi;
        const auto eps = std::abs((nn - n0) / n0);
        phi.set_eps() = eps;
        // std::cout << iteration << "_ " << eps << "\n";
        if (eps < 1.0e-7)
          break;
        n0 = nn;
      }
    }
  }

  // std::ofstream of("hf-new.txt");
  // const auto &gr = *rgrid;
  // of << "r ";
  // for (auto &psi : orbitals) {
  //   of << "\"" << psi.symbol(true) << "\" ";
  // }
  // of << "\n";
  // for (std::size_t i = 0; i < gr.num_points(); i++) {
  //   of << gr.r(i) << " ";
  //   for (auto &psi : orbitals) {
  //     of << psi.f(i) << " ";
  //   }
  //   of << "\n";
  // }

  return 0;
}

//******************************************************************************
void ContinuumOrbitals::clear() { orbitals.clear(); }
