#include "Modules/exampleModule.hpp"
#include "DiracOperator/DiracOperator.hpp" //For E1 operator
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/ostream.hpp"

namespace Module {

void exampleModule(const IO::InputBlock &, const Wavefunction &empty_wf) {
  const auto grid = empty_wf.grid();
  const auto alpha = empty_wf.alpha();
  const auto var_alpha = alpha / PhysConst::alpha;

  std::stringstream outs, outp1, outp3;

  // {1, 2, 3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88, 119, 120}
  for (const auto Z :
       {1, 3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88, 119, 120}) {

    const auto num_core_el = AtomData::z_previous_semifilled(Z);

    const auto Zion = Z - num_core_el;

    Wavefunction wf(grid.params(), Nuclear::form_nucleus(Z), var_alpha);
    std::cout << wf.nucleus() << "\n";

    wf.set_HF("HartreeFock", 0.0, AtomData::guessCoreConfigStr(num_core_el));
    wf.solve_core();
    wf.solve_valence("8sp");

    // Now, include QED:
    Wavefunction wfQ(grid.params(), Nuclear::form_nucleus(Z), var_alpha);
    wfQ.set_HF("HartreeFock", 0.0, AtomData::guessCoreConfigStr(num_core_el));
    wfQ.radiativePotential(
        IO::InputBlock{"QED", "Ueh=1.0; SE_h=0; SE_l=0; SE_m=0.0;"}, false,
        false);
    wfQ.solve_core();
    wfQ.solve_valence("8sp");

    const auto h_hf = DiracOperator::generate(
        "hfs", IO::InputBlock{"hfs", "mu=1.0; I=1.0; nuc_mag=Ball;"}, wf);

    wf.formBasis(
        SplineBasis::Parameters{"30spdfg", 40, 7, 1.0e-4, 0.0, 40.0, false});
    wfQ.formBasis(
        SplineBasis::Parameters{"30spdfg", 40, 7, 1.0e-4, 0.0, 40.0, false});

    ExternalField::TDHF dV0(h_hf.get(), wf.vHF());
    ExternalField::TDHF dVQ(h_hf.get(), wfQ.vHF());

    dV0.set_eta(0.85);
    dVQ.set_eta(0.85);

    ExternalField::DiagramRPA rpa(h_hf.get(), wf.basis(), wf.vHF(), "");
    ExternalField::DiagramRPA rpaQ(h_hf.get(), wfQ.basis(), wfQ.vHF(), "");
    if (num_core_el != 0) {
      dV0.solve_core(0.0);
      dVQ.solve_core(0.0);
      rpa.solve_core(0.0);
      rpaQ.solve_core(0.0);
    }

    for (const auto kappa : {-1, 1, -2}) {

      const auto v0 =
          *std::find_if(wf.valence().cbegin(), wf.valence().cend(),
                        [kappa](const auto &F) { return F.kappa() == kappa; });

      const auto vQ = *wfQ.getState(v0.n(), v0.kappa());
      assert(v0 == vQ);

      const auto me0 = h_hf->reducedME(v0, v0);
      const auto meQ = h_hf->reducedME(vQ, vQ);

      auto dv0 = 0.0;
      auto dvQ = 0.0;
      auto dvD0 = 0.0;
      auto dvDQ = 0.0;
      if (num_core_el != 0) {
        dv0 = dV0.dV(v0, v0);
        dvQ = dVQ.dV(vQ, vQ);
        dvD0 = rpa.dV(v0, v0);
        dvDQ = rpaQ.dV(vQ, vQ);
      }

      const auto k = 1.0 / (PhysConst::alpha / M_PI);
      const auto F_hf = k * (meQ - me0) / me0;
      const auto F_tdhf = k * (meQ + dvQ - me0 - dv0) / (me0 + dv0);
      const auto F_rpa = k * (meQ + dvDQ - me0 - dvD0) / (me0 + dvD0);

      const auto eps = std::abs((F_tdhf - F_rpa) / (F_tdhf + F_rpa));

      auto &out = kappa == -1 ? outs : kappa == 1 ? outp1 : outp3;

      fmt::print(
          out,
          "{:4s} {:3} {:4} {:11s} {:4s}  {:+.5e}  {:+.5e}  {:+.5e}  {:.0e}\n",
          wf.atomicSymbol(), Z, Zion, wf.coreConfiguration_nice(),
          v0.shortSymbol(), F_hf, F_tdhf, F_rpa, eps);
    }
  }
  std::cout << "\n\n";
  fmt::print(
      "{:4s} {:3s} {:4s} {:11s} {:4s}   {:11s}   {:11s}   {:11s}  RPA_eps\n",
      "Atom", "Z", "Zion", "Core", "val.", "F_hf", "F_tdhf", "F_diag");
  std::cout << outs.str();
  std::cout << "\n";
  std::cout << outp1.str();
  std::cout << "\n";
  std::cout << outp3.str();
  std::cout << "\n";
}

} // namespace Module
