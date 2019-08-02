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
//
#include <fstream>

namespace Module {

std::vector<double> matrixElements(const UserInputBlock &input,
                                   const Wavefunction &wf) {

  // XXX For now, only works for diagonal, just RADIAL INT!
  // XXX Need better way to deal with angular part? Output either RadInt of
  // RedMatEl

  std::string ThisModule = "MatrixElements::";
  auto operator_str = input.name().substr(ThisModule.length());

  std::cout << "\n" << ThisModule << " Operator: " << operator_str << "\n";
  auto h = generateOperator(operator_str, input, wf);

  std::vector<double> matrix_elements; // good idea?
  matrix_elements.reserve(wf.valence_orbitals.size());
  for (const auto &psi : wf.valence_orbitals) {

    auto factor = 1.0;
    if (operator_str == "hfs") { // XXX temp work-around
      double j = psi.j();
      factor = psi.k / (j * (j + 1.));
    }
    auto me = psi * (h * psi) * factor;
    matrix_elements.push_back(me);
    std::cout << psi.symbol() << " ; R = " << me << "\n";
  }
  return matrix_elements; // make this optional somehow?
} // namespace Module

//******************************************************************************
DiracOperator generateOperator(const std::string &operator_str,
                               const UserInputBlock &input,
                               const Wavefunction &wf) {
  //
  const std::string ThisModule = "MatrixElements::" + operator_str;

  if (operator_str == "hfs") {
    // XXX Lots of this (particularly F(r) part) should go into Operator!!

    auto default_mu = Nuclear::find_mu(wf.Znuc(), wf.Anuc());
    auto default_I = Nuclear::find_spin(wf.Znuc(), wf.Anuc());
    auto default_rfm = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());

    auto mu = input.get("mu", default_mu);
    auto I_nuc = input.get("I", default_I);
    auto r_rmsfm = input.get("rrms", default_rfm);
    auto r_nucfm = std::sqrt(5. / 3) * r_rmsfm;
    auto r_nucau = r_nucfm / PhysConst::aB_fm;
    auto Fr_str = input.get<std::string>("F(r)", "ball");

    std::cout << "Hyperfine structure: " << wf.atom() << "\n"
              << "Using " << Fr_str << " nuclear distro for F(r)\n"
              << "w/ mu = " << mu << ", I = " << I_nuc << ", r_N = " << r_nucfm
              << "fm = " << r_nucau << "au  (r_rms=" << r_rmsfm << "fm)\n";
    std::cout << "Points inside nucleus: " << wf.getRadialIndex(r_nucau)
              << "\n";
    std::cout << "Units: MHz\n";

    auto Fr = HyperfineOperator::sphericalBall_F();
    if (Fr_str == "shell")
      Fr = HyperfineOperator::sphericalShell_F();
    else if (Fr_str == "pointlike")
      Fr = HyperfineOperator::pointlike_F();
    else if (Fr_str == "VolotkaBW") {
      auto default_pi = Nuclear::find_parity(wf.Znuc(), wf.Anuc());
      auto pi = input.get("parity", default_pi);
      auto l_tmp = int(I_nuc + 0.5 + 0.0001);
      auto l = ((l_tmp % 2 == 0) == (pi == 1)) ? l_tmp : l_tmp - 1;
      l = input.get("l", l); // can override derived 'l' (not recommended)
      auto gl_default = wf.Znuc() % 2 == 0 ? 0 : 1; // unparied proton?
      auto gl = input.get<int>("gl", gl_default);
      std::cout << "Bohr-Weiskopf (Volotka formula) for valence";
      if (gl == 1)
        std::cout << " proton ";
      else if (gl == 0)
        std::cout << " neturon ";
      else
        std::cout << " gl=" << gl << "??? program will run, but prob wrong!\n";
      std::cout << "with l=" << l << "\n";
      Fr = HyperfineOperator::generate_F_BW(mu, I_nuc, l, gl);
    } else if (Fr_str == "doublyOddBW") {

      auto mu1 = input.get<double>("mu1");
      auto gl1 = input.get<int>("gl1"); // 1 or 0 (p or n)
      if (gl1 != 0 && gl1 != 1) {
        std::cout << "FAIL: in " << ThisModule << " " << Fr_str
                  << "; have gl1=" << gl1 << " but need 1 or 0\n";
        return DiracOperator(0);
      }
      auto l1 = input.get<double>("l1");
      auto l2 = input.get<double>("l2");
      auto I1 = input.get<double>("I1");
      auto I2 = input.get<double>("I2");

      Fr = HyperfineOperator::generate_F_BW_doublyOdd(mu, I_nuc, mu1, I1, l1,
                                                      gl1, I2, l2);
    } else if (Fr_str != "ball") {
      std::cout << "FAIL: in " << ThisModule << " " << Fr_str
                << " invalid nuclear distro. Check spelling\n";
      return DiracOperator(0);
    }

    auto unit = PhysConst::Hartree_MHz;

    auto print_FQ = input.get<bool>("printF", false);
    if (print_FQ) {
      std::ofstream of(Fr_str + ".txt");
      for (auto r : wf.rgrid.r) {
        of << r * PhysConst::aB_fm << " "
           << Fr(r * PhysConst::aB_fm, r_nucau * PhysConst::aB_fm) << "\n";
      }
    }

    // Later can add different distro's; ball is default:
    return HyperfineOperator(mu * unit, I_nuc, r_nucau, wf.rgrid, Fr);
  } else if (operator_str == "E1") {
    return E1Operator(wf.rgrid);
  } else if (operator_str == "r") {
    auto power = input.get("power", 1.0);
    std::cout << ThisModule << ": r^(" << power << ")\n";
    return RadialOperator(wf.rgrid, power);
  }

  std::cerr << "\nFAILED to find operator: " << ThisModule
            << " in generateOperator. Returning NULL operator (0)\n";
  return DiracOperator(0); // null operator
}

} // namespace Module
