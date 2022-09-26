#include "Modules/continuum.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <vector>

namespace Module {

void continuum(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"energy", "List. energy for cntm states (>0) [0.5]"},
       {"max_l", "maximum l. Will calculate orbital for each energy and each l "
                 "(and j) [0]"},
       {"filename",
        "filename to output continuum orbitals. If blank, will not write"},
       {"operator",
        "Operator to calculate matrix elements (e.g., E1) [blank by default]"},
       {"options", "options specific to operator [blank by default]"},
       {"rpa", "Include RPA (TDHF for now)? [false]"},
       {"omega", "Frequency for RPA [0.0]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto en_list = input.get("energy", std::vector{0.5});
  const auto lmax = input.get("max_l", 0);
  const auto fname = input.get("filename", std::string{""});

  // Method options for solveContinuumHF
  auto force_rescale = input.get<bool>("force_rescale", false);
  // auto subtract_self = input.get<bool>("subtract_self", false);
  auto force_orthog = input.get<bool>("force_orthog", false);

  // n=0 invalid -> nullptr
  auto n = input.get<int>("n", 0);
  auto k = input.get<int>("kappa", 0);
  // auto p_psi = wf.getState(n, k);
  auto p_psi = (n == 0) ? nullptr : wf.getState(n, k);

  if (p_psi != nullptr) {
    std::cout << "\n"
              << "State for n = " << n << " and kappa = " << k << ": "
              << p_psi->symbol() << "\n";
    // std::cout << "Bound state: " << p_psi->symbol() << "\n";
  }

  std::cout << "\nContinuum Orbitals:\n";
  auto cntm = ContinuumOrbitals(wf);

  std::cout << "For l up to l=" << lmax << "\n";
  std::cout << "At energy: ";
  for (const auto en_c : en_list) {
    std::cout << en_c << ", ";
    cntm.solveContinuumHF(en_c, lmax);
  }
  std::cout << "\n";

  //-----------------------------------------------
  std::cout << "\nCheck orthogonanilty:\n";
  cntm.check_orthog(true);

  // Matrix elements:
  //-----------------------------------------------
  const auto oper = input.get<std::string>("operator", "");
  if (oper != "") {

    // Get optional 'options' for operator
    const auto h_options = input.getBlock("options");
    const auto h = DiracOperator::generate(
        oper, h_options ? *h_options : IO::InputBlock{}, wf);

    std::cout << "\nMatrix elements of " << h->name() << "\n";

    bool eachFreqQ = false;

    const auto rpaQ = input.get("rpa", false);
    const auto omega = input.get("omega", 0.0);

    auto rpa = ExternalField::TDHF(h.get(), wf.vHF());
    const auto p_rpa = rpaQ ? &rpa : nullptr;

    const auto mes_c = ExternalField::calcMatrixElements(
        wf.core(), cntm.orbitals, h.get(), p_rpa, omega, eachFreqQ);

    const auto mes_v = ExternalField::calcMatrixElements(
        wf.valence(), cntm.orbitals, h.get(), p_rpa, omega, eachFreqQ);

    if (rpaQ) {
      std::cout << "              h(0)           h(1)           h(RPA)\n";
    } else {
      std::cout << "              h(0)\n";
    }
    std::cout << "Core:\n";
    for (const auto &me : mes_c) {
      std::cout << me << "\n";
    }
    if (!wf.valence().empty()) {
      std::cout << "Valence:\n";
    }
    for (const auto &me : mes_v) {
      std::cout << me << "\n";
    }
  }

  //-----------------------------------------------
  // Write continuum wavefunctions to file
  if (fname != "") {
    std::cout << "\nWriting orbitals to " << fname << "\n";
    std::ofstream of(fname);
    const auto &gr = wf.grid();
    of << "r ";
    for (const auto &Fc : cntm.orbitals) {
      of << "\"f" << Fc.symbol(true) << "(" << Fc.en() << ")\" ";
    }
    for (const auto &Fc : cntm.orbitals) {
      of << "\"g" << Fc.symbol(true) << "(" << Fc.en() << ")\" ";
    }
    of << "\n";
    for (std::size_t i = 0; i < gr.num_points(); i++) {
      of << gr.r(i) << " ";
      for (const auto &Fc : cntm.orbitals) {
        of << Fc.f(i) << " ";
      }
      for (const auto &Fc : cntm.orbitals) {
        of << Fc.g(i) << " ";
      }
      of << "\n";
    }
  }
}

} // namespace Module
