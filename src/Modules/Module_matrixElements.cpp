#include "Module_matrixElements.hpp"
#include "Dirac/DiracOperator.hpp"
#include "Dirac/Operators.hpp"
#include "Dirac/Wavefunction.hpp"
#include "IO/UserInput.hpp"
#include "Physics/Nuclear.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

namespace Module {

void matrixElements(const UserInputBlock &input, const Wavefunction &wf) {

  std::string ThisModule = "MatrixElements::";
  auto operator_str = input.name().substr(ThisModule.length());

  const bool radial_int = input.get("radialIntegral", false);
  auto which_str = radial_int ? "(radial integral). " : "(reduced). ";
  std::cout << "\n"
            << ThisModule << which_str << " Operator: " << operator_str << "\n";
  const auto h = generateOperator(operator_str, input, wf);

  const bool print_both = input.get("printBoth", false);
  const bool diagonal_only = input.get("onlyDiagonal", false);

  const auto units = input.get<std::string>("units", "au");
  double un = 1.0;
  if (units == "MHz")
    un = PhysConst::Hartree_MHz;

  for (const auto &phia : wf.valence_orbitals) {
    for (const auto &phib : wf.valence_orbitals) {
      if (h->isZero(phia, phib))
        continue;
      if (diagonal_only && phib != phia)
        continue;
      if (!print_both && phib < phia)
        continue;
      std::cout << h->rme_symbol(phia, phib) << ": ";
      if (radial_int)
        printf("%12.5e\n", h->radialIntegral(phia, phib) * un);
      else
        printf("%12.5e\n", h->reducedME(phia, phib));
    }
  }
} // namespace Module

//******************************************************************************
std::unique_ptr<DiracOperator> generateOperator(const std::string &operator_str,
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

    std::cout << "\nHyperfine structure: " << wf.atom() << "\n"
              << "Using " << Fr_str << " nuclear distro for F(r)\n"
              << "w/ mu = " << mu << ", I = " << I_nuc << ", r_N = " << r_nucfm
              << "fm = " << r_nucau << "au  (r_rms=" << r_rmsfm << "fm)\n";
    std::cout << "Points inside nucleus: " << wf.getRadialIndex(r_nucau)
              << "\n";

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
      std::cout << "with l=" << l << " (pi=" << pi << ")\n";
      Fr = HyperfineOperator::volotkaBW_F(mu, I_nuc, l, gl);
    } else if (Fr_str == "doublyOddBW") {

      auto mu1 = input.get<double>("mu1");
      auto gl1 = input.get<int>("gl1"); // 1 or 0 (p or n)
      if (gl1 != 0 && gl1 != 1) {
        std::cout << "FAIL: in " << ThisModule << " " << Fr_str
                  << "; have gl1=" << gl1 << " but need 1 or 0\n";
        return std::make_unique<NullOperator>(NullOperator());
      }
      auto l1 = input.get<double>("l1");
      auto l2 = input.get<double>("l2");
      auto I1 = input.get<double>("I1");
      auto I2 = input.get<double>("I2");

      Fr =
          HyperfineOperator::doublyOddBW_F(mu, I_nuc, mu1, I1, l1, gl1, I2, l2);
    } else if (Fr_str != "ball") {
      std::cout << "FAIL: in " << ThisModule << " " << Fr_str
                << " invalid nuclear distro. Check spelling\n";
      return std::make_unique<NullOperator>(NullOperator());
    }

    auto print_FQ = input.get<bool>("printF", false);
    if (print_FQ) {
      std::ofstream of(Fr_str + ".txt");
      for (auto r : wf.rgrid.r) {
        of << r * PhysConst::aB_fm << " "
           << Fr(r * PhysConst::aB_fm, r_nucau * PhysConst::aB_fm) << "\n";
      }
    }

    return std::make_unique<HyperfineOperator>(
        HyperfineOperator(mu, I_nuc, r_nucau, wf.rgrid, Fr));
  } else if (operator_str == "E1") {
    auto gauge = input.get<std::string>("gauge", "lform");
    if (gauge != "vform")
      return std::make_unique<E1Operator>(E1Operator(wf.rgrid));
    std::cout << "(v-form [velocity gauge])\n";
    return std::make_unique<E1Operator_vform>(
        E1Operator_vform(wf.rgrid, wf.get_alpha()));
  } else if (operator_str == "r") {
    auto power = input.get("power", 1.0);
    std::cout << "r^(" << power << ")\n";
    return std::make_unique<RadialFuncOperator>(
        RadialFuncOperator(wf.rgrid, power));
  } else if (operator_str == "pnc") {
    double tdflt = Nuclear::default_t; // approximate_t_skin(wf.Anuc());
    auto r_rms = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
    double cdflt = Nuclear::c_hdr_formula_rrms_t(r_rms);
    auto c = input.get("c", cdflt);
    auto t = input.get("t", tdflt);
    std::cout << "PNC [-i(Q/N)e-11]\n";
    return std::make_unique<PNCnsiOperator>(
        PNCnsiOperator(c, t, wf.rgrid, -wf.Nnuc()));
  } else if (operator_str == "M1") {
    return std::make_unique<M1Operator>(M1Operator());
  }

  std::cerr << "\nFAILED to find operator: " << ThisModule
            << " in generateOperator. Returning NULL operator (0)\n";
  return std::make_unique<NullOperator>(NullOperator());
}

} // namespace Module
