#include "ContinuumOrbitals.h"
#include "ADAMS_solveLocalBS.h"
#include "ADAMS_solveLocalContinuum.h"
#include "ATI_atomInfo.h"
#include "ElectronOrbitals.h"
#include "FPC_physicalConstants.h"
#include <cmath>
#include <string>
#include <vector>

//******************************************************************************
ContinuumOrbitals::ContinuumOrbitals(const ElectronOrbitals &wf, int izion)
/*
Initialise object:
 * Copies grid and potential info, since these must always match bound states!
If none given, will assume izion should be 1! Not 100% always!
*/
{
  // Grid:
  NGPb = wf.rgrid.ngp;
  h = wf.rgrid.du;
  r.clear();
  r = wf.rgrid.r;
  drdt.clear();
  drdt = wf.rgrid.drdu;

  alpha = wf.get_alpha();
  Z = wf.Znuc();

  // Check Zion. Will normally be 0 for neutral atom. Make -1
  double tmp_Zion =
      -1 * wf.rgrid.r[NGPb - 5] * (wf.vnuc[NGPb - 5] + wf.vdir[NGPb - 5]);

  double scale = 1;
  if (fabs(tmp_Zion - izion) > 0.01)
    scale = double(Z - izion) / (wf.rgrid.r[NGPb - 5] * wf.vdir[NGPb - 5]);

  // Local part of the potential:
  v.clear();
  v = wf.vnuc;
  if (wf.vdir.size() != 0) {
    for (int i = 0; i < NGPb; i++)
      v[i] += wf.vdir[i] * scale;
  }

  // Re-Check overal charge of atom (-1)
  // For neutral atom, should be 1 (usually, since cntm is ionisation state)
  // r->inf, v(r) = -Z_ion/r
  tmp_Zion = -1 * wf.rgrid.r[NGPb - 5] * v[NGPb - 5];
  // std::cout<<"Zion="<<tmp_Zion<<" = "<<-1*wf.rgrid.r[NGPb-5]*v[NGPb-5]<<"\n";
  Zion = izion;
}

//******************************************************************************
int ContinuumOrbitals::solveLocalContinuum(double ec, int max_l)
/*Overloaded, assumes min_l=0*/
{
  return solveLocalContinuum(ec, 0, max_l);
}

//******************************************************************************
int ContinuumOrbitals::solveLocalContinuum(double ec, int min_l, int max_l)
/*
Solved the Dirac equation for local potential for positive energy (no mc2)
continuum (un-bound) states [partial waves].
 * Goes well past NGP, looks for asymptotic region, where wf is sinosoidal
 * Uses fit to known exact H-like for normalisation.
*/
{

  // Find 'inital guess' for asymptotic region:
  double lam = 1.e7; // XXX ???
  double r_asym = (Zion + sqrt(4. * lam * ec + pow(Zion, 2))) / (2. * ec);

  // Check if 'h' is small enough for oscillating region:
  double h_target = (M_PI / 15) / sqrt(2. * ec);
  if (h > h_target) {
    std::cout << "WARNING 61 CntOrb: Grid not dense enough for ec=" << ec
              << " (h=" << h << ", need h<" << h_target << ")\n";
    if (h > 2 * h_target) {
      std::cout << "FAILURE 64 CntOrb: Grid not dense enough for ec=" << ec
                << " (h=" << h << ", need h<" << h_target << ")\n";
      return 1;
    }
  }

  // Set up temporary continuum grid:
  // Note: will have different grid sizes for different energies!
  // That's why we make a new (temporary) vector, rc
  std::vector<double> rc = r;
  std::vector<double> drdtc = drdt;
  std::vector<double> vc = v;
  int NGPc = NGPb;

  // Fill (extend) the temporary grid:
  double last_r = r[NGPb - 1];
  int i_asym = NGPb - 1;
  while (last_r < 1.2 * r_asym) {
    double r_new = last_r + h;
    if (r_new >= r_asym && last_r < r_asym)
      i_asym = NGPc - 1;
    rc.push_back(r_new);
    drdtc.push_back(1.);
    // drdtc.push_back(hc/wf.rgrid.du);
    vc.push_back(-Zion / r_new);
    NGPc++;
    last_r = r_new;
  }

  int MAX_STATES = 100;
  for (int i = 0; i < MAX_STATES; i++) { // loop through each k state
    int k = int(pow(-1, i + 1) * ceil(0.5 * (i + 1)));
    int l = (abs(2 * k + 1) - 1) / 2;
    if (l < min_l)
      continue;
    if (l > max_l)
      break;
    std::vector<double> pc(NGPb), qc(NGPb); // only as long as bound-state grid
    ADAMS::solveContinuum(pc, qc, ec, vc, k, rc, drdtc, h, NGPb, NGPc, i_asym,
                          alpha);
    f.push_back(pc);
    g.push_back(qc);
    en.push_back(ec);
    kappa.push_back(k);
  }

  return 0; // XXX code?
}

//******************************************************************************
int ContinuumOrbitals::solveZeffContinuum(double ec, double Zeff, int min_l,
                                          int max_l)
/*
Solved the Dirac equation for local potential for positive energy (no mc2)
continuum (un-bound) states [partial waves].
 * Goes well past NGP, looks for asymptotic region, where wf is sinosoidal
 * Uses fit to known exact H-like for normalisation.
XXX This no longer works! (well, it still works)
On irder to solve for Zeff, just make v use Zeff !!!!!!!
*/
{
  // Find 'inital guess' for asymptotic region:
  double lam = 1.e7;
  double r_asym = (Zeff + sqrt(4. * lam * ec + pow(Zeff, 2))) / (2. * ec);

  // Check if 'h' is small enough for oscillating region:
  double h_target = (M_PI / 15) / sqrt(2. * ec);
  if (h > h_target) {
    std::cout << "WARNING 61 CntOrb: Grid not dense enough for ec=" << ec
              << " (h=" << h << ", need h<" << h_target << ")\n";
    if (h > 2 * h_target) {
      std::cout << "FAILURE 64 CntOrb: Grid not dense enough for ec=" << ec
                << " (h=" << h << ", need h<" << h_target << ")\n";
      return 1;
    }
  }

  // Set up temporary continuum grid:
  // Note: will have different grid sizes for different energies!
  // That's why we make a new (temporary) vector, rc
  std::vector<double> rc = r;
  std::vector<double> drdtc = drdt;
  int NGPc = NGPb;

  // Fill (extend) the temporary grid:
  double last_r = r[NGPb - 1];
  int i_asym = NGPb - 1;
  while (last_r < 1.2 * r_asym) {
    double r_new = last_r + h;
    if (r_new >= r_asym && last_r < r_asym)
      i_asym = NGPc - 1;
    rc.push_back(r_new);
    drdtc.push_back(1.);
    // vc.push_back(-Zeff/r_new);
    NGPc++;
    last_r = r_new;
  }

  // Nuclear potential:
  std::vector<double> vc;
  vc.reserve(rc.size());
  for (auto rci : rc)
    vc.push_back(-Zeff / rci);

  int MAX_STATES = 100;
  for (int i = 0; i < MAX_STATES; i++) { // loop through each k state
    int k = int(pow(-1, i + 1) * ceil(0.5 * (i + 1)));
    int l = (abs(2 * k + 1) - 1) / 2;
    if (l < min_l)
      continue;
    if (l > max_l)
      break;
    std::vector<double> pc(NGPb), qc(NGPb); // only as long as bound-state grid
    ADAMS::solveContinuum(pc, qc, ec, vc, k, rc, drdtc, h, NGPb, NGPc, i_asym,
                          alpha);
    f.push_back(pc);
    g.push_back(qc);
    en.push_back(ec);
    kappa.push_back(k);
  }

  return 0;
}

//******************************************************************************
void ContinuumOrbitals::clear()
/**/
{
  f.clear();
  g.clear();
  en.clear();
  kappa.clear();
}
