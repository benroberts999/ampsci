#include "Modules/exampleModule.hpp"
#include "DiracOperator/DiracOperator.hpp" //For E1 operator
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void exampleModule(const IO::InputBlock &, const Wavefunction &wfX) {

  std::stringstream out;

  for (const auto Z : {3, 11, 19, 37, 55, 87, 119}) {

    Wavefunction wf(wfX.grid().params(), Nuclear::form_nucleus(Z));
    std::cout << wf.nucleus() << "\n";

    wf.set_HF("HartreeFock", 0.0, AtomData::guessCoreConfigStr(Z - 1));
    wf.solve_core();
    wf.solve_valence("8sp");

    // Now, include QED:
    Wavefunction wfQ(wfX.grid().params(), Nuclear::form_nucleus(Z));
    wfQ.set_HF("HartreeFock", 0.0, AtomData::guessCoreConfigStr(Z - 1));
    wfQ.radiativePotential(
        IO::InputBlock{"QED", "Ueh=1.0; SE_h=0; SE_l=0; SE_m=0.0;"}, false,
        false);
    wfQ.solve_core();
    wfQ.solve_valence("8sp");

    const auto h_hf = DiracOperator::generate(
        "hfs", IO::InputBlock{"hfs", "mu=1.0; I=1.0; nuc_mag=pointlike;"}, wf);

    const auto v0 = wf.valence().front();
    const auto vQ = wfQ.valence().front();

    ExternalField::TDHF dV0(h_hf.get(), wf.vHF());
    dV0.solve_core(0.0);
    ExternalField::TDHF dVQ(h_hf.get(), wfQ.vHF());
    dVQ.solve_core(0.0);

    const auto me0 = h_hf->reducedME(v0, v0);
    const auto meQ = h_hf->reducedME(vQ, vQ);

    const auto dv0 = dV0.dV(v0, v0);
    const auto dvQ = dVQ.dV(vQ, vQ);

    out << Z << " " << (meQ - me0) / me0 << " "
        << (meQ + dvQ - me0 - dv0) / (me0 + dv0) << "\n";
  }
  std::cout << out.str();
}

} // namespace Module
