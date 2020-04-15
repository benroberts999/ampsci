#include "Modules/Module_matrixElements.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "DiracOperator/Operators.hpp"
#include "HF/ExternalField.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/RadiativePotential.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>

namespace Module {

//******************************************************************************
void matrixElements(const IO::UserInputBlock &input, const Wavefunction &wf) {

  std::string ThisModule = "MatrixElements::";
  auto operator_str = input.name().substr(ThisModule.length());
  // nb: "check" is done in 'generate operator'

  const bool radial_int = input.get("radialIntegral", false);

  // spacial case: HFS A (MHz)
  bool AhfsQ = (operator_str == "hfs" && !radial_int);

  auto which_str =
      radial_int ? "(radial integral). " : AhfsQ ? " A (MHz). " : "(reduced). ";

  auto h = generateOperator(operator_str, input, wf);
  std::cout << "\n"
            << ThisModule << which_str << " Operator: " << h->name() << "\n";
  std::cout << "Units: " << h->units() << "\n";

  const bool print_both = input.get("printBoth", false);
  const bool diagonal_only = input.get("onlyDiagonal", false);

  const bool rpaQ = input.get("rpa", true);

  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = str_om == "each" || str_om == "Each";
  const auto omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  const auto factor = input.get("factor", 1.0);
  if (factor != 1.0)
    std::cout << "With extra factor: " << factor << "\n";

  // XXX Always same kappa for get_Vlocal?
  auto rpa = HF::ExternalField(h.get(), wf.core, wf.get_Vlocal(), wf.alpha);
  std::unique_ptr<HF::ExternalField> rpa0; // for first-order

  if (h->freqDependantQ && !eachFreqQ)
    h->updateFrequency(omega);

  if (!eachFreqQ && rpaQ) {
    rpa.solve_TDHFcore(omega, 1, false);
    rpa0 =
        std::make_unique<HF::ExternalField>(rpa); // store first-order snapshot
    rpa.solve_TDHFcore(omega);
  } else {
    rpa0 = std::make_unique<HF::ExternalField>(rpa); // Solved later
  }

  // Fb -> Fa = <a||h||b>
  for (const auto &Fb : wf.valence) {
    for (const auto &Fa : wf.valence) {

      if (h->isZero(Fa.k, Fb.k))
        continue;
      if (diagonal_only && Fb != Fa)
        continue;
      if (!print_both && Fb > Fa)
        continue;
      const auto ww = std::abs(Fa.en - Fb.en);
      if (eachFreqQ && h->freqDependantQ) {
        h->updateFrequency(ww);
      }
      if (eachFreqQ && rpaQ) {
        rpa0->clear_dPsi();
        rpa0->solve_TDHFcore(ww, 1, false); // wastes a little time
        rpa.clear_dPsi();                   // in case last one didn't work!
        rpa.solve_TDHFcore(ww); // re-solve at new frequency (not from scratch)
      }
      std::cout << h->rme_symbol(Fa, Fb) << ": ";
      // Special case: HFS A:
      auto a = AhfsQ ? DiracOperator::Hyperfine::convertRMEtoA(Fa, Fb) : 1.0;
      a *= factor;
      if (radial_int) {
        printf("%13.6e\n", h->radialIntegral(Fa, Fb));
      } else if (rpaQ) {
        auto dV = rpa.dV_ab(Fa, Fb);
        auto dV0 = rpa0->dV_ab(Fa, Fb);
        printf("%13.6e  %13.6e  %13.6e\n", h->reducedME(Fa, Fb) * a,
               (h->reducedME(Fa, Fb) + dV0) * a,
               (h->reducedME(Fa, Fb) + dV) * a);
      } else {
        printf("%13.6e\n", h->reducedME(Fa, Fb) * a);
      }
    }
  }
}

//******************************************************************************
void calculateBohrWeisskopf(const IO::UserInputBlock &input,
                            const Wavefunction &wf)
//
{
  using namespace DiracOperator;

  IO::UserInputBlock point_in("MatrixElements::hfs", input);
  IO::UserInputBlock ball_in("MatrixElements::hfs", input);
  IO::UserInputBlock BW_in("MatrixElements::hfs", input);
  point_in.add("F(r)=pointlike");
  ball_in.add("F(r)=ball");
  if (wf.Anuc() % 2 == 0)
    BW_in.add("F(r)=doublyOddBW");
  else
    BW_in.add("F(r)=VolotkaBW");

  auto hp = generateOperator("hfs", point_in, wf);
  auto hb = generateOperator("hfs", ball_in, wf);
  auto hw = generateOperator("hfs", BW_in, wf);

  std::cout << "\nTabulate A (Mhz), and Bohr-Weisskopf effect eps(%): "
            << wf.atom()
            << "\n       |A:      Point         Ball           SP |e:    "
               "Ball         SP\n";
  for (const auto &phi : wf.valence) {
    auto Ap = Hyperfine::hfsA(hp.get(), phi);
    auto Ab = Hyperfine::hfsA(hb.get(), phi);
    auto Aw = Hyperfine::hfsA(hw.get(), phi);
    auto Fball = ((Ab / Ap) - 1.0) * 100.0; //* M_PI * PhysConst::c;
    auto Fbw = ((Aw / Ap) - 1.0) * 100.0;   //* M_PI * PhysConst::c;
    printf("%7s| %12.5e %12.5e %12.5e | %9.5f  %9.5f \n", phi.symbol().c_str(),
           Ap, Ab, Aw, Fball, Fbw);
  }
}

//******************************************************************************
std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(const std::string &operator_str,
                 const IO::UserInputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  //
  const std::string ThisModule = "MatrixElements::" + operator_str;

  std::vector<std::string> check_list = {
      "radialIntegral", "printBoth", "onlyDiagonal", "units", "rpa",
      "omega",          "factor"};
  auto jointCheck = [&](const std::vector<std::string> &in) {
    check_list.insert(check_list.end(), in.begin(), in.end());
    return check_list;
  };

  if (operator_str == "hfs") {
    input.checkBlock(
        jointCheck({"mu", "I", "rrms", "F(r)", "parity", "l", "gl", "mu1",
                    "gl1", "l1", "l2", "I1", "I2", "printF"}));

    auto isotope = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
    auto mu = input.get("mu", isotope.mu);
    auto I_nuc = input.get("I", isotope.I_N);
    auto r_rmsfm = input.get("rrms", isotope.r_rms);
    auto r_nucfm = std::sqrt(5. / 3) * r_rmsfm;
    auto r_nucau = r_nucfm / PhysConst::aB_fm;
    auto Fr_str = input.get<std::string>("F(r)", "ball");

    std::cout << "\nHyperfine structure: " << wf.atom() << "\n"
              << "Using " << Fr_str << " nuclear distro for F(r)\n"
              << "w/ mu = " << mu << ", I = " << I_nuc << ", r_N = " << r_nucfm
              << "fm = " << r_nucau << "au  (r_rms=" << r_rmsfm << "fm)\n";
    std::cout << "Points inside nucleus: " << wf.rgrid.getIndex(r_nucau)
              << "\n";

    auto Fr = Hyperfine::sphericalBall_F();
    if (Fr_str == "shell")
      Fr = Hyperfine::sphericalShell_F();
    else if (Fr_str == "pointlike")
      Fr = Hyperfine::pointlike_F();
    else if (Fr_str == "VolotkaBW") {
      auto pi = input.get("parity", isotope.parity);
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
      Fr = Hyperfine::volotkaBW_F(mu, I_nuc, l, gl);
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

      Fr = Hyperfine::doublyOddBW_F(mu, I_nuc, mu1, I1, l1, gl1, I2, l2);
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

    return std::make_unique<Hyperfine>(
        Hyperfine(mu, I_nuc, r_nucau, wf.rgrid, Fr));
  } else if (operator_str == "E1") {
    input.checkBlock(jointCheck({"gauge"}));
    auto gauge = input.get<std::string>("gauge", "lform");
    if (gauge != "vform")
      return std::make_unique<E1>(wf.rgrid);
    // std::cout << "(v-form [velocity gauge])\n";
    return std::make_unique<E1_vform>(wf.alpha, 0.0);
  } else if (operator_str == "Ek") {
    input.checkBlock(jointCheck({"k"}));
    auto k = input.get("k", 1);
    return std::make_unique<Ek>(Ek(wf.rgrid, k));
  } else if (operator_str == "r") {
    input.checkBlock(jointCheck({"power"}));
    auto power = input.get("power", 1.0);
    std::cout << "r^(" << power << ")\n";
    return std::make_unique<RadialF>(RadialF(wf.rgrid, power));
  } else if (operator_str == "pnc") {
    input.checkBlock(jointCheck({"c", "t"}));
    double tdflt = Nuclear::default_t; // approximate_t_skin(wf.Anuc());
    auto r_rms = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
    double cdflt = Nuclear::c_hdr_formula_rrms_t(r_rms);
    auto c = input.get("c", cdflt);
    auto t = input.get("t", tdflt);
    std::cout << "PNC [-i(Q/N)e-11]\n";
    return std::make_unique<PNCnsi>(PNCnsi(c, t, wf.rgrid, -wf.Nnuc()));
  } else if (operator_str == "M1") {
    return std::make_unique<M1>(wf.rgrid, wf.alpha, 0.0);
  } else if (operator_str == "Hrad_el") {
    input.checkBlock(
        jointCheck({"Simple", "Ueh", "SE_h", "SE_l", "rcut", "scale_rN"}));
    const auto x_Simple = input.get("Simple", 0.0);
    const auto x_Ueh = input.get("Ueh", 1.0);
    const auto x_SEe_h = input.get("SE_h", 1.0);
    const auto x_SEe_l = input.get("SE_l", 1.0);
    // const auto x_SEm = input.get("SE_m", 1.0);
    const auto rcut = input.get("rcut", 5.0);
    const auto scale_rN = input.get("scale_rN", 1.0);
    const auto r_rms_Fermi = scale_rN * wf.get_nuclearParameters().r_rms;
    const auto Hel = RadiativePotential::form_Hel(wf.rgrid.r, x_Simple, x_Ueh,
                                                  x_SEe_h, x_SEe_l, r_rms_Fermi,
                                                  wf.Znuc(), wf.alpha, rcut);
    return std::make_unique<Hrad_el>(Hrad_el(Hel));
  } else if (operator_str == "Hrad_mag") {
    input.checkBlock(jointCheck({"SE_m", "rcut", "scale_rN"}));
    const auto x_SEm = input.get("SE_m", 1.0);
    const auto rcut = input.get("rcut", 5.0);
    const auto scale_rN = input.get("scale_rN", 1.0);
    const auto r_rms_Fermi = scale_rN * wf.get_nuclearParameters().r_rms;
    const auto Hmag = RadiativePotential::form_Hmag(
        wf.rgrid.r, x_SEm, r_rms_Fermi, wf.Znuc(), wf.alpha, rcut);
    return std::make_unique<Hrad_mag>(Hrad_mag(Hmag));
  }

  std::cerr << "\nFAILED to find operator: " << ThisModule
            << " in generateOperator. Returning NULL operator (0)\n";
  return std::make_unique<NullOperator>(NullOperator());
}

//******************************************************************************
void calculateLifetimes(const IO::UserInputBlock &input,
                        const Wavefunction &wf) {
  std::cout << "\nLifetimes:\n";

  input.checkBlock({"E1", "E2"});
  auto doE1 = input.get("E1", true);
  auto doE2 = input.get("E2", false);
  if (doE1 && !doE2)
    std::cout << "Including E1 only.\n";
  if (!doE1 && doE2)
    std::cout << "Including E2 only.\n";

  DiracOperator::E1 he1(wf.rgrid);
  DiracOperator::Ek he2(wf.rgrid, 2);
  auto alpha = wf.alpha;
  auto alpha3 = alpha * alpha * alpha;
  auto alpha2 = alpha * alpha;
  auto dVE1 = HF::ExternalField(&he1, wf.core, wf.get_Vlocal(), alpha);
  auto dVE2 = HF::ExternalField(&he2, wf.core, wf.get_Vlocal(), alpha);

  auto to_s = PhysConst::time_s;

  for (const auto &Fa : wf.valence) {
    std::cout << "\n" << Fa.symbol() << "\n";
    auto Gamma = 0.0;

    if (doE1) {
      for (const auto &Fn : wf.valence) {
        if (Fn.en >= Fa.en || he1.isZero(Fn.k, Fa.k))
          continue;
        auto w = Fa.en - Fn.en;
        dVE1.solve_TDHFcore(w, 40);
        auto d = he1.reducedME(Fn, Fa) + dVE1.dV_ab(Fn, Fa);
        auto g_n = (4.0 / 3) * w * w * w * d * d / (Fa.twojp1());
        Gamma += g_n;
        std::cout << "  E1 --> " << Fn.symbol() << ": ";
        printf("w=%7.5f, |d|=%7.5f, g=%10.4eau\n", w, std::abs(d),
               g_n * alpha3);
      }
    }
    if (doE2) {
      for (const auto &Fn : wf.valence) {
        if (Fn.en >= Fa.en || he2.isZero(Fn.k, Fa.k))
          continue;
        auto w = Fa.en - Fn.en;
        dVE2.solve_TDHFcore(w, 40);
        auto d = he2.reducedME(Fn, Fa) + dVE2.dV_ab(Fn, Fa);
        auto g_n = (1.0 / 15) * w * w * w * w * w * d * d / (Fa.twojp1());
        Gamma += g_n * alpha2;
        std::cout << "  E2 --> " << Fn.symbol() << ": ";
        printf("w=%7.5f, |q|=%7.5f, g=%10.4eau\n", w, std::abs(d),
               g_n * alpha3);
      }
    }

    printf("Gamma = %10.4eau = %10.4e/s\n", Gamma * alpha3,
           Gamma * alpha3 / to_s);
    printf("tau = %10.4es\n", to_s / alpha3 / Gamma);
  }
}

} // namespace Module
