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

  for (const auto Z :
       {1,  2,  3,  4,  5,  6,  11, 12, 13, 14, 19, 20,  21,  22,  37,
        38, 39, 40, 55, 56, 57, 58, 87, 88, 89, 90, 119, 120, 121, 122}) {

    const auto num_core_el = Z == 2 ? 0 : AtomData::z_previous_semifilled(Z);

    const auto Zion = Z - num_core_el;

    auto nucleus = Nuclear::form_nucleus(Z);
    Wavefunction wf(grid.params(), nucleus, var_alpha);

    std::cout << "\n------------------------"
                 "\n------------------------"
                 "\n------------------------\n\n";
    std::cout << wf.atom() << "\n";
    std::cout << wf.nucleus() << "\n\n";

    wf.set_HF("HartreeFock", 0.0, AtomData::guessCoreConfigStr(num_core_el));
    wf.solve_core();

    const auto n_val = DiracSpinor::max_n(wf.core()) + 2;
    const auto val = std::to_string(n_val) + "sp";
    std::cout << "Val " << val << "\n";

    wf.solve_valence(val);

    // Now, include QED:
    Wavefunction wfQ(grid.params(), nucleus, var_alpha);
    wfQ.set_HF("HartreeFock", 0.0, AtomData::guessCoreConfigStr(num_core_el));
    wfQ.radiativePotential(
        IO::InputBlock{"QED", "Ueh=1.0; SE_h=0; SE_l=0; SE_m=0.0;"}, false,
        true);
    wfQ.solve_core();
    wfQ.solve_valence(val);

    wf.printCore();
    std::cout << "No QED:\n";
    wf.printValence();
    std::cout << "Vac Pol:\n";
    wfQ.printValence();

    const auto h_hf = DiracOperator::generate(
        "hfs", IO::InputBlock{"hfs", "mu=1.0; I=1.0; nuc_mag=pointlike;"}, wf);

    wf.formBasis(SplineBasis::Parameters{"40spdfg", 45, 9, 1.0e-4, 0.0,
                                         45.0 / std::sqrt(Zion), false});
    wfQ.formBasis(SplineBasis::Parameters{"40spdfg", 45, 9, 1.0e-4, 0.0,
                                          45.0 / std::sqrt(Zion), false});

    ExternalField::TDHF dV0(h_hf.get(), wf.vHF());
    ExternalField::TDHF dVQ(h_hf.get(), wfQ.vHF());

    dV0.set_eta(0.90);
    dVQ.set_eta(0.90);

    ExternalField::DiagramRPA rpa(h_hf.get(), wf.basis(), wf.vHF(), "");
    ExternalField::DiagramRPA rpaQ(h_hf.get(), wfQ.basis(), wfQ.vHF(), "");
    if (num_core_el != 0) {
      dV0.solve_core(0.0, 1);
      dVQ.solve_core(0.0, 1);
      rpa.solve_core(0.0, 1);
      rpaQ.solve_core(0.0, 1);
    }

    for (const auto kappa : {-1, 1, -2}) {

      const auto pv0 =
          std::find_if(wf.valence().cbegin(), wf.valence().cend(),
                       [kappa](const auto &F) { return F.kappa() == kappa; });
      if (pv0 == wf.valence().cend()) {
        std::cout << __LINE__ << "\n";
        continue;
      }
      auto v0 = *pv0;

      const auto vQ = *wfQ.getState(v0.n(), v0.kappa());
      assert(v0 == vQ);

      const auto me0 = h_hf->reducedME(v0, v0);
      const auto meQ = h_hf->reducedME(vQ, vQ);

      const auto A = DiracOperator::Hyperfine::convert_RME_to_AB(1, v0.kappa(),
                                                                 v0.kappa());

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

      fmt::print(out,
                 "{:4s} {:3} {:4} {:4s} {:4s} {:.3f}  {:+.5e}  {:+.5e}  "
                 "{:+.5e}  {:+.5e}  {:+.5e} "
                 " {:.0e}\n",
                 wf.atomicSymbol(), Z, Zion, wf.coreConfiguration_nice(),
                 v0.shortSymbol(), wf.get_rrms(), A * me0, A * (me0 + dvD0),
                 F_hf, F_tdhf, F_rpa, eps);
    }
  }
  std::cout << "\n\n";
  fmt::print("{:4s} {:3s} {:4s} {:4s} {:4s} {:5s}   {:11s}   {:11s}   {:11s}   "
             "{:11s}   {:11s}  "
             "RPA_eps\n",
             "Atom", "Z", "Zion", "Core", "val.", "Rrms", "A_hf", "A_rpa",
             "F_hf", "F_tdhf", "F_diag");
  std::cout << outs.str();
  std::cout << "\n";
  std::cout << outp1.str();
  std::cout << "\n";
  std::cout << outp3.str();
  std::cout << "\n";
}

} // namespace Module
