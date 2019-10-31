#include "Dirac/ContinuumOrbitals.hpp"
#include "Adams/DiracODE.hpp"
#include "Dirac/Wavefunction.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomInfo.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <cmath>
#include <string>
#include <vector>

//******************************************************************************
ContinuumOrbitals::ContinuumOrbitals(const Wavefunction &wf, int izion)
    : p_rgrid(&wf.rgrid), Z(wf.Znuc()), Zion(izion), alpha(wf.get_alpha())
// Initialise object:
//  * Copies grid and potential info, since these must always match bound
//  states!
// If none given, will assume izion should be 1! Not 100% always!
{
  auto num_pointsb = wf.rgrid.num_points;

  // Check Zion. Will normally be 0 for neutral atom. Make -1
  double tmp_Zion = -1 * wf.rgrid.r[num_pointsb - 5] *
                    (wf.vnuc[num_pointsb - 5] + wf.vdir[num_pointsb - 5]);

  // Note: Because I don't include exchange, need to re-scale the potential
  // Assumes z_ion=1 (i.e., one electron ejected from otherwise neutral atom)
  // This is equivilent to using the averaged hartree potential
  // NOTE: I _could_ have izion = (tmp_Zion +1). But easier for now to just
  // input z_ion (sometimes, might want to do something different)
  double scale = 1;
  if (std::fabs(tmp_Zion - izion) > 0.01)
    scale = double(Z - izion) /
            (wf.rgrid.r[num_pointsb - 5] * wf.vdir[num_pointsb - 5]);

  // Local part of the potential:
  v.clear();
  v = wf.vnuc;
  if (wf.vdir.size() != 0) {
    for (auto i = 0ul; i < num_pointsb; i++) {
      v[i] += wf.vdir[i] * scale;
    }
  }

  // Re-Check overal charge of atom (-1)
  // For neutral atom, should be 1 (usually, since cntm is ionisation state)
  // r->inf, v(r) = -Z_ion/r
  tmp_Zion = -1 * wf.rgrid.r[num_pointsb - 5] * v[num_pointsb - 5];

  if (std::fabs(tmp_Zion - Zion) > 0.01) {
    std::cout << "\nWARNING: [cntm] Zion incorrect?? Is this OK??\n";
    std::cout << "Zion=" << tmp_Zion << " = "
              << -1 * wf.rgrid.r[num_pointsb - 5] * v[num_pointsb - 5] << " "
              << izion << "\n";
    std::cin.get();
  }
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
  double lam = 1.0e7; // XXX ???
  double r_asym =
      (Zion + std::sqrt(4. * lam * ec + std::pow(Zion, 2))) / (2. * ec);

  // Check if 'h' is small enough for oscillating region:
  double h_target = (M_PI / 15) / std::sqrt(2. * ec);
  auto h = p_rgrid->du;
  if (h > h_target) {
    std::cout << "WARNING 61 CntOrb: Grid not dense enough for ec=" << ec
              << " (h=" << h << ", need h<" << h_target << ")\n";
    if (h > 2 * h_target) {
      std::cout << "FAILURE 64 CntOrb: Grid not dense enough for ec=" << ec
                << " (h=" << h << ", need h<" << h_target << ")\n";
      return 1;
    }
  }

  ExtendedGrid cgrid(*p_rgrid, 1.2 * r_asym);

  auto vc = v;
  for (auto i = p_rgrid->num_points; i < cgrid.num_points; i++) {
    vc.push_back(-Zion / cgrid.r[i]);
  }

  for (int i = 0; true; ++i) { // loop through each k state
    auto k = AtomInfo::kappaFromIndex(i);
    auto l = AtomInfo::l_k(k);
    if (l < min_l)
      continue;
    if (l > max_l)
      break;

    // guess as asymptotic region:
    auto i_asym = cgrid.getIndex(r_asym); // - 1;

    DiracSpinor phi(0, k, *p_rgrid);
    phi.en = ec;
    DiracODE::solveContinuum(phi, vc, cgrid, i_asym, alpha);

    orbitals.push_back(phi);
  }

  return 0;
}

//******************************************************************************
void ContinuumOrbitals::clear() { orbitals.clear(); }
