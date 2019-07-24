#include "Module_matrixElements.hpp"
#include "DiracOperator.hpp"
#include "HartreeFockClass.hpp"
#include "Operators.hpp"
#include "Physics/Nuclear.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "UserInput.hpp"
#include "Wavefunction.hpp"
#include <iostream>
#include <string>

namespace Module {

void matrixElements(const std::string &module, const UserInput &input,
                    const Wavefunction &wf) {

  // XXX For now, only works for diagonal, just RADIAL INT!

  std::string ThisModule = "MatrixElements::";
  auto operator_str = module.substr(ThisModule.length());

  auto h = generateOperator(operator_str, input, wf);

  for (const auto &psi : wf.valence_orbitals) {
    // XXX temp: ///
    double j = psi.j();
    auto factor = psi.k / (j * (j + 1.));
    std::cout << psi.symbol() << " ; R = " << psi * (h * psi) * factor << "\n";
  }
}

//******************************************************************************
DiracOperator generateOperator(const std::string &operator_str,
                               const UserInput &input, const Wavefunction &wf) {
  //
  static const std::string ThisModule = "MatrixElements::" + operator_str;

  if (operator_str == "hfs") {
    // XXX Lots of this (particularly F(r) part) should go into Operator!!

    auto default_mu = Nuclear::find_mu(wf.Znuc(), wf.Anuc());
    auto default_I = Nuclear::find_spin(wf.Znuc(), wf.Anuc());
    auto default_rfm = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
    // auto default_r = wf.get_rrms();

    auto mu = input.get(ThisModule, "mu", default_mu);
    auto I_nuc = input.get(ThisModule, "I", default_I);
    auto r_nucfm = input.get(ThisModule, "r_rms", default_rfm);
    auto r_nucau = r_nucfm / PhysConst::aB_fm;
    auto Fr_str = input.get<std::string>(ThisModule, "F(r)", "ball");

    std::cout << "\nHyperfine structure:\n"
              << "Using " << Fr_str << " nuclear distro for F(r)\n"
              << "w/ mu = " << mu << ", I = " << I_nuc << " rrms = " << r_nucfm
              << "fm = " << r_nucau << "au\n";
    std::cout << "Units: MHz\n";

    auto Fr = HyperfineOperator::sphericalBall_F();
    if (Fr_str == "shell")
      Fr = HyperfineOperator::sphericalShell_F();
    if (Fr_str == "pointlike")
      Fr = HyperfineOperator::pointlike_F();
    if (Fr_str == "VolotkaBW") {
      auto l = input.get<double>(ThisModule, "l");
      auto gl = input.get<int>(ThisModule, "gl");
      std::cout << "Bohr-Weiskopf (Volotka formula) for valence";
      if (gl == 1)
        std::cout << " proton ";
      else if (gl == 0)
        std::cout << " neturon ";
      else
        std::cout << " gl=" << gl << "??? program will run, but prob wrong!\n";
      std::cout << "with l=" << l << "\n";
      Fr = HyperfineOperator::generate_F_BW(mu, I_nuc, l, double(gl));
    }
    auto unit = PhysConst::Hartree_MHz;

    // Later can add different distro's; ball is default:
    return HyperfineOperator(mu * unit, I_nuc, r_nucau, wf.rgrid, Fr);
  } else if (operator_str == "E1") {
    return E1Operator(wf.rgrid);
  } else if (operator_str == "r") {
    auto power = input.get(ThisModule, "power", 1.0);
    return RadialOperator(wf.rgrid, power);
  }

  std::cerr << "\nFAILED to find operator: " << ThisModule
            << " in generateOperator. Returning NULL operator (0)\n";
  return DiracOperator(0.0); // null operator
}

} // namespace Module
