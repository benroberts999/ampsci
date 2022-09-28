#include "Modules/qed.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/matrixElements.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/String.hpp"
#include <numeric>
#include <string>
#include <vector>

namespace Module {

void QED(const IO::InputBlock &input, const Wavefunction &wf) {

  // Check input options for spelling mistakes etc.:
  input.check(
      {{"rcut", "Maximum radius (au) to calculate Rad Pot for [5.0]"},
       {"scale_rN", "Scale factor for Nuclear size. 0 for pointlike, 1 for "
                    "typical [1.0]"},
       {"scale_l", "List of doubles. Extra scaling factor for each l e.g., "
                   "1,0,1 => include for s and d, but not for p [1.0]"},
       {"core_QED", "Inlcude QED into core (or only valence)? [true]"},
       {"MatrixElements{}", "For "
                            "QED corrects to MEs. Input block; takes mostly "
                            "same inputs an Module::MatrixElements."}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // We should not run this module if outer wavefunction aleady includes QED, as this will double-count!
  if (wf.vrad() != nullptr) {
    std::cout << "\nDo not run QED module if QED corrections are inlcuded at "
                 "main level!\n"
              << "Module not run\n\n";
    return;
  }

  // We simply 'copy' the correlation potential across from the
  // outer-wavefunction. This is usually fine, but means QED was not included
  // in the basis/Green's function that was used to construct Sigma.
  if (wf.Sigma() != nullptr) {
    std::cout
        << "\nNote: QED corrections not included into correlations\n"
        << "This is probably fine, but should be confirmed independently\n";
  }

  // This module calculates QED corections:
  // 1. Energy corrections, without relaxation
  // 2. Energy corrections, with relaxation
  // 3. Optional: corrections to matrix elements

  // Read in options for QED:
  const auto rcut = input.get("rcut", 5.0);
  const auto scale_rn = input.get("scale_rN", 1.0);
  const auto r_N_au =
      std::sqrt(5.0 / 3.0) * wf.nucleus().r_rms() * scale_rn / PhysConst::aB_fm;
  const auto x_spd = input.get("scale_l", std::vector{1.0});
  bool include_qed_core = input.get("core_QED", true);

  // Construct the radiative potentials (seperately, and a 'total' one)
  // nb: don't write to disk, allows easier update of, e.g., rcut
  std::cout << "\n";
  auto Ueh = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                         {1.0, 0.0, 0.0, 0.0, 0.0}, x_spd, true, false);
  auto SEh = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                         {0.0, 1.0, 0.0, 0.0, 0.0}, x_spd, true, false);
  auto SEl = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                         {0.0, 0.0, 1.0, 0.0, 0.0}, x_spd, true, false);
  auto SEm = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                         {0.0, 0.0, 0.0, 1.0, 0.0}, x_spd, true, false);
  auto WK = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                        {0.0, 0.0, 0.0, 0.0, 1.0}, x_spd, true, false);
  auto Vtot = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                          {1.0, 1.0, 1.0, 1.0, 1.0}, x_spd, false, false);

  //----------------------------------------------------------------------------

  // Create operators, used to find first-order energy corrections
  DiracOperator::Vrad Vu(Ueh);
  DiracOperator::Vrad Vh(SEh);
  DiracOperator::Vrad Vl(SEl);
  DiracOperator::Vrad Vm(SEm);
  DiracOperator::Vrad Vw(WK);
  DiracOperator::Vrad Vt(Vtot);

  // Include QED into Hartree-Fock:
  // We 'copy' everything over from outer wavefunction, then add QED potential.
  Wavefunction wf_u = wf;
  Wavefunction wf_h = wf;
  Wavefunction wf_l = wf;
  Wavefunction wf_m = wf;
  Wavefunction wf_w = wf;
  Wavefunction wf_t = wf;
  if (wf.vHF()) {
    wf_u.vHF()->set_Vrad(Ueh);
    wf_h.vHF()->set_Vrad(SEh);
    wf_l.vHF()->set_Vrad(SEl);
    wf_m.vHF()->set_Vrad(SEm);
    wf_w.vHF()->set_Vrad(WK);
    wf_t.vHF()->set_Vrad(Vtot);
  } else {
    std::cout << "?";
    return;
  }

  // If including QED into core, we will then re-solve HF for the core.
  if (include_qed_core) {
    std::cout << "\nRe-solving Hartree-Fock for core states, including:\n";
    std::cout << "Uehling potential:\n";
    wf_u.solve_core();
    std::cout << "Self-energy (high frequency):\n";
    wf_h.solve_core();
    std::cout << "Self-energy (low frequency):\n";
    wf_l.solve_core();
    std::cout << "Self-energy (magnetic):\n";
    wf_m.solve_core();
    std::cout << "Wichmann-Kroll (approximate):\n";
    wf_w.solve_core();
    std::cout << "Total radiative potential:\n";
    wf_t.solve_core();

    std::cout << "\nCore energies (including total QED):\n";
    wf_t.printCore(false);
  } else {
    std::cout
        << "\nIncluding QED radiative potential into valence states only\n";
  }

  // We clear the valence wavefunctions, and re-solve them from scratch.
  wf_u.valence().clear();
  wf_h.valence().clear();
  wf_l.valence().clear();
  wf_m.valence().clear();
  wf_w.valence().clear();
  wf_t.valence().clear();

  // Solve for valence states with QED:
  // nb: hartreeFockBrueckner() has no effect if Sigma doesn't exist
  const auto valence_list = DiracSpinor::state_config(wf.valence());
  std::cout << "\nRe-solving Hartree-Fock for valence states " << valence_list
            << ", inclduding:\n";

  std::cout << "Uehling potential:\n";
  wf_u.solve_valence(valence_list);
  wf_u.hartreeFockBrueckner();

  std::cout << "Self-energy (high frequency):\n";
  wf_h.solve_valence(valence_list);
  wf_h.hartreeFockBrueckner();

  std::cout << "Self-energy (low frequency):\n";
  wf_l.solve_valence(valence_list);
  wf_l.hartreeFockBrueckner();

  std::cout << "Self-energy (magnetic):\n";
  wf_m.solve_valence(valence_list);
  wf_m.hartreeFockBrueckner();

  std::cout << "Wichmann-Kroll (approximate):\n";
  wf_w.solve_valence(valence_list);
  wf_w.hartreeFockBrueckner();

  std::cout << "Total radiative potential:\n";
  wf_t.solve_valence(valence_list);
  wf_t.hartreeFockBrueckner();

  std::cout << "\nValence energies (including total QED):\n";
  wf_t.printValence(false);

  //----------------------------------------------------------------------------
  // Output:

  // 1. First-order QED corrections to energies
  std::cout << "\nFirst-order QED corrections to energies (/cm):\n";
  std::cout
      << "State  Uehl      SE(h)     SE(l)     SE(m)     WK        Total\n";
  for (auto &v : wf.valence()) {
    const auto deu = Vu.radialIntegral(v, v) * PhysConst::Hartree_invcm;
    const auto deh = Vh.radialIntegral(v, v) * PhysConst::Hartree_invcm;
    const auto del = Vl.radialIntegral(v, v) * PhysConst::Hartree_invcm;
    const auto dem = Vm.radialIntegral(v, v) * PhysConst::Hartree_invcm;
    const auto dew = Vw.radialIntegral(v, v) * PhysConst::Hartree_invcm;
    const auto det = Vt.radialIntegral(v, v) * PhysConst::Hartree_invcm;
    // const auto det2 = deu + deh + del + dem + dew;
    printf("%4s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n", v.shortSymbol().c_str(),
           deu, deh, del, dem, dew, det);
  }

  // 2. Total QED corrections to energies
  std::cout << "\nTotal QED corrections to energies (/cm), with relaxation:\n";
  std::cout
      << "State  Uehl      SE(h)     SE(l)     SE(m)     WK        Total\n";
  for (auto &v0 : wf.valence()) {
    const auto &vu = *wf_u.getState(v0.n(), v0.kappa());
    const auto &vh = *wf_h.getState(v0.n(), v0.kappa());
    const auto &vl = *wf_l.getState(v0.n(), v0.kappa());
    const auto &vm = *wf_m.getState(v0.n(), v0.kappa());
    const auto &vw = *wf_w.getState(v0.n(), v0.kappa());
    const auto &vt = *wf_t.getState(v0.n(), v0.kappa());

    const auto deu = (vu.en() - v0.en()) * PhysConst::Hartree_invcm;
    const auto deh = (vh.en() - v0.en()) * PhysConst::Hartree_invcm;
    const auto del = (vl.en() - v0.en()) * PhysConst::Hartree_invcm;
    const auto dem = (vm.en() - v0.en()) * PhysConst::Hartree_invcm;
    const auto dew = (vw.en() - v0.en()) * PhysConst::Hartree_invcm;
    const auto det = (vt.en() - v0.en()) * PhysConst::Hartree_invcm;
    printf("%4s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
           v0.shortSymbol().c_str(), deu, deh, del, dem, dew, det);
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Matrix elements:

  const auto me_input = input.getBlock("MatrixElements");
  // If we aren't given this block, we're done:
  if (me_input == std::nullopt)
    return;

  std::cout << "\n";
  IO::print_line();
  std::cout << "QED corrections to matrix elements\n";

  me_input->check({{"operator", "e.g., E1, hfs"},
                   {"options{}", "options specific to operator; blank by dflt"},
                   {"rpa", "true or false (only TDHF method) [true]"},
                   {"omega", "Frequency for RPA [0.0]"}});

  const auto oper = me_input->get<std::string>("operator", "");
  // Get optional 'options' for operator
  // auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = me_input->getBlock("options");
  // if (tmp_opt) {
  //   h_options = *tmp_opt;
  // }
  const auto h =
      DiracOperator::generate(oper, tmp_opt ? *tmp_opt : IO::InputBlock{}, wf);

  const auto rpaQ = me_input->get("rpa", true);
  const auto omega = me_input->get("omega", 0.0);

  if (h->parity() == 1 && rpaQ) {
    std::cout << "\n\n*CAUTION*:\n RPA (TDHF method) may not work for this "
                 "operator.\n Consider using diagram or basis method\n\n";
  }

  // set up RPA, for each QED sub-operator:
  using namespace ExternalField;
  std::unique_ptr<CorePolarisation> dV0{nullptr};
  std::unique_ptr<CorePolarisation> dVu{nullptr};
  std::unique_ptr<CorePolarisation> dVh{nullptr};
  std::unique_ptr<CorePolarisation> dVl{nullptr};
  std::unique_ptr<CorePolarisation> dVm{nullptr};
  std::unique_ptr<CorePolarisation> dVw{nullptr};
  std::unique_ptr<CorePolarisation> dVt{nullptr};
  if (rpaQ) {
    std::cout << "Including RPA: ";
    std::cout << "TDHF method\n";
    dV0 = std::make_unique<TDHF>(h.get(), wf.vHF());
    dVu = std::make_unique<TDHF>(h.get(), wf_u.vHF());
    dVh = std::make_unique<TDHF>(h.get(), wf_h.vHF());
    dVl = std::make_unique<TDHF>(h.get(), wf_l.vHF());
    dVm = std::make_unique<TDHF>(h.get(), wf_m.vHF());
    dVw = std::make_unique<TDHF>(h.get(), wf_w.vHF());
    dVt = std::make_unique<TDHF>(h.get(), wf_t.vHF());
  }

  std::cout << "\nCalculating matrix elements ";
  if (rpaQ) {
    std::cout << "(and solving RPA) ";
  }
  std::cout << "including:\n";
  std::cout << "No QED:\n";
  const auto me0 = calcMatrixElements(wf.valence(), h.get(), dV0.get(), omega);

  std::cout << "Uehling:\n";
  const auto meu =
      calcMatrixElements(wf_u.valence(), h.get(), dVu.get(), omega);

  std::cout << "Self-energy (high-frequency):\n";
  const auto meh =
      calcMatrixElements(wf_h.valence(), h.get(), dVh.get(), omega);

  std::cout << "Self-energy (low-frequency):\n";
  const auto mel =
      calcMatrixElements(wf_l.valence(), h.get(), dVl.get(), omega);

  std::cout << "Self-energy (magnetic):\n";
  const auto mem =
      calcMatrixElements(wf_m.valence(), h.get(), dVm.get(), omega);

  std::cout << "Wichmann-Kroll:\n";
  const auto mew =
      calcMatrixElements(wf_w.valence(), h.get(), dVw.get(), omega);

  std::cout << "Total radiative potential:\n";
  const auto met =
      calcMatrixElements(wf_t.valence(), h.get(), dVt.get(), omega);

  // Print matrix elements, without QED:
  std::cout << "\nReduced matrix elements (" << h->name() << "), no QED:\n";
  if (rpaQ) {
    std::cout << "             h(0)           h(1)           h(RPA)\n";
  }
  for (const auto &me : me0) {
    std::cout << me << "\n";
  }

  // Print matrix elements, with total QED:
  std::cout << "\nReduced matrix elements (" << h->name()
            << "), including full QED:\n";
  if (rpaQ) {
    std::cout << "             h(0)           h(1)           h(RPA)\n";
  }
  for (const auto &me : met) {
    std::cout << me << "\n";
  }

  // QED corrections to MEs (no RPA):
  std::cout << "\nQED corrections to reduced matrix elements (no RPA), in au\n";
  std::cout << "States     Uehl       SE(h)      SE(l)      SE(m)      WK     "
               "    Total\n";
  for (std::size_t i = 0; i < me0.size(); ++i) {
    const auto [a, b, hab, d1, dv] = me0.at(i);
    printf("%4s %4s ", a.c_str(), b.c_str());
    for (const auto p : {&meu, &meh, &mel, &mem, &mew, &met}) {
      const auto [qa, qb, qh, qd1, qdv] = (*p).at(i);
      assert(qa == a && qb == b);
      auto del0 = qh - hab;
      // auto deldv = (qh + qdv) - (hab + dv);
      printf("%10.3e ", del0);
    }
    std::cout << "\n";
  }

  // QED corrections to MEs (with RPA):
  if (rpaQ) {
    // QED corrections to MEs:
    std::cout << "\nQED corrections to matrix elements (with RPA), in au\n";
    std::cout
        << "States     Uehl       SE(h)      SE(l)      SE(m)      WK     "
           "    Total\n";
    for (std::size_t i = 0; i < me0.size(); ++i) {
      const auto [a, b, hab, d1, dv] = me0.at(i);
      printf("%4s %4s ", a.c_str(), b.c_str());
      for (const auto p : {&meu, &meh, &mel, &mem, &mew, &met}) {
        const auto [qa, qb, qh, qd1, qdv] = (*p).at(i);
        assert(qa == a && qb == b);
        // auto del0 = qh - hab;
        auto deldv = (qh + qdv) - (hab + dv);
        printf("%10.3e ", deldv);
      }
      std::cout << "\n";
    }
  }
}

} // namespace Module
