#include "HF/Breit.hpp"
#include "DiracOperator/DiracOperator.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "Modules/exampleModule.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include <string>

namespace Module {

void Breit(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check({{"", "Currently takes no options"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto Br = HF::Breit{};

  // To separate Gaus/Retardation parts
  auto G = HF::Breit{};
  auto R = HF::Breit{};
  G.update_scale(1.0, 1.0, 1.0, 0.0, 0.0);
  R.update_scale(1.0, 0.0, 0.0, 1.0, 1.0);

  for (const auto &v : wf.valence()) {
    std::cout << v << " " << v * Br.VbrFa(v, wf.core()) << "\n";
  }

  const auto eFermi = wf.FermiLevel();
  const auto holes =
      qip::select_if(wf.basis(), [=](auto &a) { return a.en() < eFermi; });
  const auto excited =
      qip::select_if(wf.basis(), [=](auto &a) { return a.en() > eFermi; });

  std::cout << "\n";

  std::cout << "Core energy corrections (without relaxation):\n";
  std::cout << "      E(HF)       Gaunt    Ret.    Total  \n";
  for (const auto &a : wf.core()) {
    auto e0 = a.en();
    auto eg = a * G.VbrFa(a, wf.core());
    auto er = a * R.VbrFa(a, wf.core());
    printf("%4s %10.4f %7.4f %7.4f  %7.4f\n", a.shortSymbol().c_str(), e0, eg,
           er, eg + er);
  }

  if (!wf.valence().empty()) {
    std::cout << "\nValence energy corrections (without relaxation):\n";
    std::cout << "      E(HF)     de(1)      de(2,Z1)   de(2,Z2)   de(2)      "
                 "Total\n";
  }
  for (const auto &v : wf.valence()) {
    auto e0 = v.en();
    auto de0 = v * Br.VbrFa(v, wf.core());
    auto deHF = Br.de2_HF(v, holes, excited);
    auto de2 = Br.de2(v, holes, excited);
    // std::cout << v << " " << de0 << " " << deHF + de2 << "\n";
    printf("%4s %9.6f %10.3e %10.3e %10.3e %10.3e %10.3e\n",
           v.shortSymbol().c_str(), e0, de0, deHF, de2, deHF + de2,
           de0 + deHF + de2);
  }

  std::cout << "\n";
  std::cout << "Solve Hartree-Fock again, including Breit\n";
  auto wf2 = wf;
  wf2.solve_core("HartreeFock", 1.0, wf.coreConfiguration());
  wf2.valence().clear();
  wf2.solve_valence(DiracSpinor::state_config(wf.valence()));
  wf2.printValence();

  if (!wf.valence().empty())
    std::cout << "\nValence energy corrections, including relaxation (HF "
                 "level)\n";
  for (const auto &v : wf.valence()) {
    auto v2 = wf2.getState(v.symbol());
    std::cout << v << " " << v2->en() - v.en() << "\n";
  }
}

} // namespace Module
