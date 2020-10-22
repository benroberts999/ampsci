#include "Modules/matrixElements.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "DiracOperator/Operators.hpp"
#include "HF/ExternalField.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/DiagramRPA.hpp"
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
#include <vector>

namespace Module {

//******************************************************************************
void matrixElements(const IO::UserInputBlock &input, const Wavefunction &wf) {

  input.checkBlock({"operator", "options", "rpa", "rpa_diagram", "omega",
                    "radialIntegral", "printBoth", "onlyDiagonal", "units",
                    "A_vertex", "b_vertex"});

  const auto oper = input.get<std::string>("operator", "");
  const auto h_options =
      IO::UserInputBlock(oper, input.get<std::string>("options", ""));
  const auto h = generateOperator(oper, h_options, wf, true);

  const bool radial_int = input.get("radialIntegral", false);

  // spacial case: HFS A (MHz)
  const bool AhfsQ = (oper == "hfs" && !radial_int);

  const auto which_str =
      radial_int ? " (radial integral)." : AhfsQ ? " A (MHz)." : " (reduced).";

  std::cout << "\n"
            << input.name() << which_str << " Operator: " << h->name() << "\n";
  std::cout << "Units: " << h->units() << "\n";

  const bool print_both = input.get("printBoth", false);
  const bool diagonal_only = input.get("onlyDiagonal", false);

  const bool rpaDQ = input.get("rpa_diagram", false);
  const bool rpaQ = input.get("rpa", !rpaDQ);

  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = str_om == "each" || str_om == "Each";
  const auto omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  if (rpaQ && (h->parity() == 1)) {
    std::cout << "\n\n*CAUTION*:\n RPA (TDHF method) may not work for this "
                 "operator.\n Consider using rpa_diagram method\n\n";
  }

  // XXX Make this its own module! XXX
  // Vertex QED term:
  std::unique_ptr<DiracOperator::VertexQED> hVertexQED = nullptr;
  const auto A_vertex = input.get("A_vertex", 0.0);
  if (A_vertex != 0.0) {
    const auto b_vertex = input.get("b_vertex", 1.0);
    std::cout << "Including effective vertex QED with: A=" << A_vertex
              << ", b=" << b_vertex << "\n";
    hVertexQED = std::make_unique<DiracOperator::VertexQED>(h.get(), *wf.rgrid,
                                                            A_vertex, b_vertex);
  }

  // XXX Always same kappa for get_Vlocal?
  auto rpa = HF::ExternalField(h.get(), wf.getHF());
  std::unique_ptr<HF::ExternalField> rpa0; // for first-order
  std::unique_ptr<MBPT::DiagramRPA> rpaD;

  if (h->freqDependantQ && !eachFreqQ)
    h->updateFrequency(omega);

  if (rpaQ) {
    if (!eachFreqQ) {
      rpa.solve_TDHFcore(omega, 1, false);
      // store first-order snapshot:
      rpa0 = std::make_unique<HF::ExternalField>(rpa);
      rpa.solve_TDHFcore(omega);
    } else {
      rpa0 = std::make_unique<HF::ExternalField>(rpa); // Solved later
    }
  }
  if (rpaDQ) {
    rpaD = std::make_unique<MBPT::DiagramRPA>(h.get(), wf.basis, wf.core,
                                              wf.identity());
    if (!eachFreqQ)
      rpaD->rpa_core(omega);
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
        if (rpa.get_eps() > 1.0e-5)
          rpa.clear_dPsi();     // in case last one didn't work!
        rpa.solve_TDHFcore(ww); // re-solve at new frequency
      }
      if (eachFreqQ && rpaDQ) {
        if (rpaD->get_eps() > 1.0e-5)
          rpaD->clear_tam(); // in case last one didn't work!
        rpaD->rpa_core(ww);
      }

      // Special case: HFS A:
      const auto a = AhfsQ ? DiracOperator::Hyperfine::convertRMEtoA(Fa, Fb)
                           : radial_int ? 1.0 / h->angularF(Fa.k, Fb.k) : 1.0;

      const auto symb =
          radial_int ? h->R_symbol(Fa, Fb) + " " : h->rme_symbol(Fa, Fb);

      if (!rpaQ && !rpaDQ) {
        std::cout << symb << ": ";
        printf("%13.6e \n", h->reducedME(Fa, Fb) * a);
      }
      if (rpaQ) {
        std::cout << symb << ": ";
        printf("%13.6e ", h->reducedME(Fa, Fb) * a);
        auto dV = rpa.dV(Fa, Fb);
        auto dV0 = rpa0->dV(Fa, Fb);
        printf(" %13.6e  %13.6e\n", (h->reducedME(Fa, Fb) + dV0) * a,
               (h->reducedME(Fa, Fb) + dV) * a);
      }
      if (rpaDQ) {
        std::cout << symb << "D:";
        printf("%13.6e ", h->reducedME(Fa, Fb) * a);
        auto dV = rpaD->dV(Fa, Fb);
        auto dV0 = rpaD->dV(Fa, Fb, true);
        printf(" %13.6e  %13.6e\n", (h->reducedME(Fa, Fb) + dV0) * a,
               (h->reducedME(Fa, Fb) + dV) * a);
      }
      // XXX Make its own module!
      if (hVertexQED) { //
        printf("   QED vertex: ");
        printf("%13.6e \n", hVertexQED->reducedME(Fa, Fb) * a);
        // Add RPA? Might be important?
      }
    }
  }
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

  DiracOperator::E1 he1(*(wf.rgrid));
  DiracOperator::Ek he2(*(wf.rgrid), 2);
  auto alpha = wf.alpha;
  auto alpha3 = alpha * alpha * alpha;
  auto alpha2 = alpha * alpha;
  auto dVE1 = HF::ExternalField(&he1, wf.getHF());
  auto dVE2 = HF::ExternalField(&he2, wf.getHF());

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
        auto d = he1.reducedME(Fn, Fa) + dVE1.dV(Fn, Fa);
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
        auto d = he2.reducedME(Fn, Fa) + dVE2.dV(Fn, Fa);
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

//******************************************************************************
//******************************************************************************
std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(const IO::UserInputBlock &input, const Wavefunction &wf,
                 bool print) {
  using namespace DiracOperator;

  for (const auto &[name, generator] : operator_list) {
    // (void)generator;
    if (/*"MatrixElements::" + name == input.name() ||*/ name == input.name())
      return generator(input, wf, print);
  }

  std::cerr << "\nFAILED to find operator: " << input.name()
            << " in generateOperator.\n";

  std::cout << "Available operators:\n";
  for (const auto &[name, generator] : operator_list) {
    (void)generator;
    std::cout << "  " << name << "\n";
  }
  std::cout << "\n";
  std::cout << "Returning NULL operator (0)\n";

  return std::make_unique<NullOperator>(NullOperator());
}

std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(std::string_view oper_name, const IO::UserInputBlock &input,
                 const Wavefunction &wf, bool print) {
  using namespace DiracOperator;

  for (const auto &[name, generator] : operator_list) {
    // (void)generator;
    if (name == oper_name)
      return generator(input, wf, print);
  }

  std::cerr << "\nFAILED to find operator: " << input.name()
            << " in generateOperator.\n";

  std::cout << "Available operators:\n";
  for (const auto &[name, generator] : operator_list) {
    (void)generator;
    std::cout << "  " << name << "\n";
  }
  std::cout << "\n";
  std::cout << "Returning NULL operator (0)\n";

  return std::make_unique<NullOperator>(NullOperator());
}

//******************************************************************************
//******************************************************************************

// inline auto jointCheck(const std::vector<std::string> &in) {
//   std::vector<std::string> check_list = {
//       "radialIntegral", "printBoth", "onlyDiagonal", "units",   "rpa",
//       "rpa_diagram",    "omega",     "A_vertex",     "b_vertex"};
//   check_list.insert(check_list.end(), in.begin(), in.end());
//   return check_list;
// }

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_E1(const IO::UserInputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock({"gauge"});
  auto gauge = input.get<std::string>("gauge", "lform");
  if (gauge != "vform")
    return std::make_unique<E1>(*(wf.rgrid));
  // std::cout << "(v-form [velocity gauge])\n";
  return std::make_unique<E1_vform>(wf.alpha, 0.0);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_Ek(const IO::UserInputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock({"k"});
  auto k = input.get("k", 1);
  return std::make_unique<Ek>(*(wf.rgrid), k);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_M1(const IO::UserInputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock({});
  return std::make_unique<M1>(*(wf.rgrid), wf.alpha, 0.0);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_hfs(const IO::UserInputBlock &input, const Wavefunction &wf,
             bool print) {
  using namespace DiracOperator;
  input.checkBlock({"mu", "I", "rrms", "F(r)", "parity", "l", "gl", "mu1",
                    "gl1", "l1", "l2", "I1", "I2", "printF", "screening"});
  auto isotope = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  auto mu = input.get("mu", isotope.mu);
  auto I_nuc = input.get("I", isotope.I_N);
  auto r_rmsfm = input.get("rrms", isotope.r_rms);
  auto r_nucfm = std::sqrt(5. / 3) * r_rmsfm;
  auto r_nucau = r_nucfm / PhysConst::aB_fm;
  auto Fr_str = input.get<std::string>("F(r)", "ball");

  if (print) {
    std::cout << "\nHyperfine structure: " << wf.atom() << "\n"
              << "Using " << Fr_str << " nuclear distro for F(r)\n"
              << "w/ mu = " << mu << ", I = " << I_nuc << ", r_N = " << r_nucfm
              << "fm = " << r_nucau << "au  (r_rms=" << r_rmsfm << "fm)\n";
    std::cout << "Points inside nucleus: " << wf.rgrid->getIndex(r_nucau)
              << "\n";
  }

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
    if (print) {
      std::cout << "Bohr-Weiskopf (Volotka formula) for valence";
      if (gl == 1)
        std::cout << " proton ";
      else if (gl == 0)
        std::cout << " neturon ";
      else
        std::cout << " gl=" << gl << "??? program will run, but prob wrong!\n";
      std::cout << "with l=" << l << " (pi=" << pi << ")\n";
    }
    Fr = Hyperfine::volotkaBW_F(mu, I_nuc, l, gl);
  } else if (Fr_str == "doublyOddBW") {

    auto mu1 = input.get<double>("mu1");
    auto gl1 = input.get<int>("gl1"); // 1 or 0 (p or n)
    if (gl1 != 0 && gl1 != 1) {
      std::cout << "FAIL: in " << input.name() << " " << Fr_str
                << "; have gl1=" << gl1 << " but need 1 or 0\n";
      return std::make_unique<NullOperator>(NullOperator());
    }
    auto l1 = input.get<double>("l1");
    auto l2 = input.get<double>("l2");
    auto I1 = input.get<double>("I1");
    auto I2 = input.get<double>("I2");

    Fr = Hyperfine::doublyOddBW_F(mu, I_nuc, mu1, I1, l1, gl1, I2, l2);
  } else if (Fr_str != "ball") {
    std::cout << "FAIL: in " << input.name() << " " << Fr_str
              << " invalid nuclear distro. Check spelling\n";
    return std::make_unique<NullOperator>();
  }

  auto print_FQ = input.get<bool>("printF", false);
  if (print_FQ) {
    std::ofstream of(Fr_str + ".txt");
    for (auto r : wf.rgrid->r) {
      of << r * PhysConst::aB_fm << " "
         << Fr(r * PhysConst::aB_fm, r_nucau * PhysConst::aB_fm) << "\n";
    }
  }

  return std::make_unique<Hyperfine>(
      Hyperfine(mu, I_nuc, r_nucau, *(wf.rgrid), Fr));
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_r(const IO::UserInputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock({"power"});
  auto power = input.get("power", 1.0);
  std::cout << "r^(" << power << ")\n";
  return std::make_unique<RadialF>(*(wf.rgrid), power);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_pnc(const IO::UserInputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock({"c", "t"});
  const auto r_rms = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
  const auto c = input.get("c", Nuclear::c_hdr_formula_rrms_t(r_rms));
  const auto t = input.get("t", Nuclear::default_t);
  return std::make_unique<PNCnsi>(c, t, *(wf.rgrid), 1.0, "iQwe-11");
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_Hrad_el(const IO::UserInputBlock &input, const Wavefunction &wf,
                 bool) {
  using namespace DiracOperator;
  input.checkBlock({"Simple", "Ueh", "SE_h", "SE_l", "rcut", "scale_rN"});
  const auto x_Simple = input.get("Simple", 0.0);
  const auto x_Ueh = input.get("Ueh", 1.0);
  const auto x_SEe_h = input.get("SE_h", 1.0);
  const auto x_SEe_l = input.get("SE_l", 1.0);
  // const auto x_SEm = input.get("SE_m", 1.0);
  const auto rcut = input.get("rcut", 5.0);
  const auto scale_rN = input.get("scale_rN", 1.0);
  const auto r_rms_Fermi = scale_rN * wf.get_nuclearParameters().r_rms;
  const auto Hel = RadiativePotential::form_Hel(wf.rgrid->r, x_Simple, x_Ueh,
                                                x_SEe_h, x_SEe_l, r_rms_Fermi,
                                                wf.Znuc(), wf.alpha, rcut);
  return std::make_unique<Hrad_el>(Hel);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_Hrad_mag(const IO::UserInputBlock &input, const Wavefunction &wf,
                  bool) {
  using namespace DiracOperator;
  input.checkBlock({"SE_m", "rcut", "scale_rN"});
  const auto x_SEm = input.get("SE_m", 1.0);
  const auto rcut = input.get("rcut", 50.0);
  const auto scale_rN = input.get("scale_rN", 1.0);
  const auto r_rms_Fermi = scale_rN * wf.get_nuclearParameters().r_rms;
  const auto Hmag = RadiativePotential::form_Hmag(
      wf.rgrid->r, x_SEm, r_rms_Fermi, wf.Znuc(), wf.alpha, rcut);
  return std::make_unique<Hrad_mag>(Hmag);
}

//******************************************************************************

} // namespace Module
