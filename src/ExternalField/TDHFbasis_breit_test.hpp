#pragma once
#include "DiracOperator/DiracOperator.hpp"
#include "TDHFbasis.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

namespace UnitTest {

//******************************************************************************
namespace helper {
std::pair<double, std::string>
do_dV_breit_basis(const Wavefunction &wf, const Wavefunction &wfB,
                  const DiracOperator::TensorOperator *const h,
                  double omega = 0.0, int max_l = 99);
} // namespace helper

//******************************************************************************
//! Compared Breit contribution to dV using TDHF and TDHFbasis methods
bool TDHFbasis_breit(std::ostream &obuff) { //
  bool pass = true;

  // Test of internal consistancy, not validation of Breit correction

  // nb: don't use many grid points, since we don't need overall accuracy.
  // Need fair splines though, since splineTDHF compared to TDHF

  // Compares the Breit correction to dV using TDHF and TDHF basis methods.
  // NB: also a test of BASIS+Breit

  std::string valence = "7s6p5d";
  SplineBasis::Parameters bspl_param;
  {
    bspl_param.states = "40spd30f";
    bspl_param.n = 50;
    bspl_param.k = 7;
    bspl_param.r0 = 1.0e-5;
    bspl_param.reps = 1.0e-6;
    bspl_param.rmax = 40.0;
  }

  // Create wavefunction object, solve HF for core+valence
  Wavefunction wf({1000, 1.0e-6, 120.0, 50.0, "loglinear", -1.0},
                  {"Cs", 133, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Xe]"); // no breit
  wf.solve_valence(valence);
  wf.formBasis(bspl_param);
  // wf.printValence();

  // Again, including Breit:
  Wavefunction wfB(wf.rgrid->params(), wf.get_nuclearParameters(), 1.0);
  wfB.solve_core("HartreeFock", 1.0, "[Xe]"); // w/Breit
  wfB.solve_valence(valence);
  wfB.formBasis(bspl_param);

  const auto hE1 = DiracOperator::E1(*wf.rgrid);
  const auto hPNC =
      DiracOperator::PNCnsi(Nuclear::c_hdr_formula_rrms_t(wf.get_rrms()),
                            Nuclear::default_t, *wf.rgrid);

  using namespace helper;
  const auto [eE1, sE1] = do_dV_breit_basis(wf, wfB, &hE1, 0.0);
  const auto [eE1w, sE1w] = do_dV_breit_basis(wf, wfB, &hE1, 0.05);
  const auto max_l = 1; // don't include d-states for pnc?
  const auto [ePNC, sPNC] = do_dV_breit_basis(wf, wfB, &hPNC, 0.0, max_l);

  pass &=
      qip::check_value(&obuff, "Breit TDHF(bs) E1 w=0 " + sE1, eE1, 0.0, 0.005);
  pass &= qip::check_value(&obuff, "Breit TDHF(bs) E1 w=0.05 " + sE1w, eE1w,
                           0.0, 0.05);
  pass &=
      qip::check_value(&obuff, "Breit TDHF(bs) PNC " + sPNC, ePNC, 0.0, 0.1);

  return pass;
}

//******************************************************************************
std::pair<double, std::string>
helper::do_dV_breit_basis(const Wavefunction &wf, const Wavefunction &wfB,
                          const DiracOperator::TensorOperator *const h,
                          double omega, int max_l) {

  auto dV_tdhf = ExternalField::TDHF(h, wf.getHF());
  auto dV_basis = ExternalField::TDHFbasis(h, wf.getHF(), wf.basis);
  // w/ Breit:
  auto dVB_tdhf = ExternalField::TDHF(h, wfB.getHF());
  auto dVB_basis = ExternalField::TDHFbasis(h, wfB.getHF(), wfB.basis);

  // double omega = 0.0;
  const auto max_iterations = 50;
  const auto print_details = true;
  std::cout << "No Breit:\n";
  dV_tdhf.solve_core(omega, max_iterations, print_details);
  std::cout << "With Breit:\n";
  dVB_tdhf.solve_core(omega, max_iterations, print_details);
  std::cout << "No Breit:\n";
  dV_basis.solve_core(omega, max_iterations, print_details);
  std::cout << "With Breit:\n";
  dVB_basis.solve_core(omega, max_iterations, print_details);

  auto eps = 0.0;
  std::string worst = "";

  std::cout << "dBr(%)  TDHF       TDHF(basis) " << h->name() << "\n";
  for (const auto &Fv : wf.valence) {
    for (const auto &Fw : wf.valence) {
      if (Fw > Fv || h->isZero(Fv.k, Fw.k) || Fv.l() > max_l || Fw.l() > max_l)
        continue;

      // lookup states in Breit wf (do this way, so no depend on order)
      const auto &FvB = *wfB.getState(Fv.n, Fv.k);
      const auto &FwB = *wfB.getState(Fw.n, Fw.k);

      const auto dv = dV_tdhf.dV(Fv, Fw);
      const auto dvB = dVB_tdhf.dV(FvB, FwB);
      const auto dv_basis = dV_basis.dV(Fv, Fw);
      const auto dvB_basis = dVB_basis.dV(FvB, FwB);

      const auto dBr = (dvB - dv) / dv;
      const auto dBr_basis = (dvB_basis - dv_basis) / dv_basis;
      const auto str = Fv.shortSymbol() + "-" + Fw.shortSymbol();
      printf("%7s %9.4e %9.4e\n", str.c_str(), dBr * 100.0, dBr_basis * 100.0);
      const auto eps_this = std::abs((dBr - dBr_basis) / (dBr + dBr_basis));
      if (eps_this > eps) {
        eps = eps_this;
        worst = str;
      }
    }
  }
  return {eps, worst};
}

} // namespace UnitTest
