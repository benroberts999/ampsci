#include "Modules/CI_Pol.hpp"
#include "CI/CI_Integrals.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/include.hpp" //For E1 operator
#include "ExternalField/TDHF.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void CI_Pol(const IO::InputBlock &input, const Wavefunction &wf) {
  // This module computes matrix elements and energies for two valence systems

  //these quantities are then used to compute second order amplitudes using the sum over states (SOS) method

  // In this example, we will solve a new wavefunction, assuming a different
  // nuclear charge distribution, and see the difference in the energies and E1
  // matrix elements this produces.

  // Read in some optional input options: A, rms, and type
  // "checkBlock" is optional, but helps prevent errors on input options:
  input.check(
      {{"J", "Angular momentum of wavefunction"},
       {"parity", "Parity of state"},
       {"operator_1", "e.g., E1, hfs (see ampsci -o for available operators)"},
       {"state_number", "specify which state with angular momentum J and "
                        "parity PI is required"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  //external field operator------------------------------------------------------------
  //get ranks of operators
  const auto oper1 = input.get<std::string>("operator_1", "");

  // Get optional 'options' for operator
  auto h1_options = IO::InputBlock(oper1, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h1_options = *tmp_opt;
  }
  const auto h1 = DiracOperator::generate(oper1, h1_options, wf);
  const auto k1 = h1->rank();
  const auto p1 = h1->parity();
  //---------------------------------------------------------------------------------------
  //states
  const auto J = input.get("J", 0);
  const auto parity = input.get("parity", 1);
  const auto nv = input.get("state_mumber", 0ul);
  //get CI wavefunction for "initial state"
  const auto wfV = wf.CIwf(J, parity);

  std::vector<DiracSpinor> orbitals = wf.basis();
  //store and fill table of single orbtial matrix elements
  ExternalField::TDHF tdhf(h1.get(), wf.vHF());
  tdhf.solve_core(0.0);
  Coulomb::meTable<double> sTable =
      ExternalField::me_table(orbitals, h1.get(), &tdhf);

  //Number of states for final allowed angular momentum and parity obtained when solving CI+MBPT

  // get allowed states for matrix elements of oprator of rank K

  //for polarisability
  double pol = 0.0;
  double pol_factor = 2.0 / (3.0 * (2.0 * J + 1.0));

  for (int J2 = J - k1; J2 <= J + k1; J2++) {
    if (J2 < 0) {
      continue;
    }
    if (J2 == 0 && J == 0) {
      continue;
    }
    const auto wfn = wf.CIwf(J2, -parity);
    const auto num = wfn->num_solutions();
    double conj_phase = Angular::neg1pow_2(2 * J - 2 * J2);
    for (std::size_t i = 0; i < num; i++) {
      double CI_ME_1 = CI::ReducedME(*wfV, nv, *wfn, i, sTable, k1, p1);
      double CI_ME_2 = CI::ReducedME(*wfn, i, *wfV, nv, sTable, k1, p1);
      double deltaE = wfn->energy(i) - wfV->energy(nv);

      std::cout << i << " " << parity << " " << CI_ME_1 << "  " << CI_ME_2
                << " " << deltaE << "\n";
      double pol_cont = CI_ME_1 * conj_phase * CI_ME_2 / deltaE;
      pol += pol_cont;
    }
  }
  //wfV.info(nv).config;
  const auto cf = (*wfV).info(nv).config;
  pol = pol_factor * pol;
  std::cout << "Polarisability of state with configuration " << cf
            << ", angular momentum " << J << " and parity " << parity << " is "
            << pol << " a_0^3\n";
}
} // namespace Module
