#include "Modules/matrixElements.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/StructureRad.hpp"
#include "Maths/Interpolator.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/String.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace Module {

//******************************************************************************
void matrixElements(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check({{"operator", "e.g., E1, hfs"},
               {"options", "options specific to operator; blank by dflt"},
               {"rpa", "true(=TDHF), false, TDHF, basis, diagram"},
               {"omega", "freq. for RPA"},
               {"radialIntegral", "false by dflt (means red. ME)"},
               {"printBoth", "print <a|h|b> and <b|h|a> (dflt false)"},
               {"onlyDiagonal", "only <a|h|a> (dflt false)"}});

  const auto oper = input.get<std::string>("operator", "");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

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

  const auto rpa_method_str = input.get("rpa", std::string("TDHF"));
  auto rpa_method = (rpa_method_str == "TDHF" || rpa_method_str == "true") ?
                        ExternalField::method::TDHF :
                        (rpa_method_str == "basis") ?
                        ExternalField::method::basis :
                        (rpa_method_str == "diagram") ?
                        ExternalField::method::diagram :
                        ExternalField::method::none;
  if (wf.core.empty())
    rpa_method = ExternalField::method::none;
  const auto rpaQ = rpa_method != ExternalField::method::none;
  // const auto rpaDQ = rpa_method == ExternalField::method::diagram;

  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = str_om == "each" || str_om == "Each";
  const auto omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  if (h->freqDependantQ) {
    std::cout << "Frequency-dependent operator; at omega = ";
    if (eachFreqQ)
      std::cout << "each transition frequency\n";
    else
      std::cout << omega << "\n";
  }

  if ((h->parity() == 1) && rpa_method == ExternalField::method::TDHF) {
    std::cout << "\n\n*CAUTION*:\n RPA (TDHF method) may not work for this "
                 "operator.\n Consider using diagram or basis method\n\n";
  }

  std::unique_ptr<ExternalField::CorePolarisation> rpa{nullptr};
  if (rpaQ)
    std::cout << "Including RPA: ";
  if (rpa_method == ExternalField::method::TDHF) {
    std::cout << "TDHF method\n";
    rpa = std::make_unique<ExternalField::TDHF>(h.get(), wf.getHF());
  } else if (rpa_method == ExternalField::method::basis) {
    std::cout << "TDHF/basis method\n";
    rpa = std::make_unique<ExternalField::TDHFbasis>(h.get(), wf.getHF(),
                                                     wf.basis);
  } else if (rpa_method == ExternalField::method::diagram) {
    std::cout << "diagram method\n";
    rpa = std::make_unique<ExternalField::DiagramRPA>(h.get(), wf.basis,
                                                      wf.core, wf.identity());
  }

  const auto mes = ExternalField::calcMatrixElements(
      wf.valence, h.get(), rpa.get(), omega, eachFreqQ, diagonal_only,
      print_both, radial_int);

  if (rpaQ) {
    std::cout << "              h(0)           h(1)           h(RPA)\n";
  }
  for (const auto &me : mes) {
    std::cout << me << "\n";
  }
}

//****************************************************************************
// Calculates Structure Radiation + Normalisation of States
void structureRad(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check({{"operator", "e.g., E1, hfs"},
               {"options", "options specific to operator; blank by dflt"},
               {"rpa", "true(=TDHF), false, TDHF, basis, diagram"},
               {"omega", "freq. for RPA"},
               {"printBoth", "print <a|h|b> and <b|h|a> (dflt false)"},
               {"onlyDiagonal", "only <a|h|a> (dflt false)"},
               {"n_minmax", "list; min,max n for core/excited: (1,inf)dflt"},
               {"splineLegs", "Use splines for diagram legs (false dflt)"}});

  // Get input options:
  const auto oper = input.get<std::string>("operator", "E1");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  const auto h = generateOperator(oper, h_options, wf, true);

  // Use spline states as diagram legs?
  const auto spline_legs = input.get("splineLegs", false);

  const auto print_both = input.get("printBoth", false);
  const auto only_diagonal = input.get("onlyDiagonal", false);
  const auto rpaQ = input.get("rpa", true);
  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = str_om == "each" || str_om == "Each";
  const auto const_omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  // min/max n (for core/excited basis)
  const auto n_minmax = input.get("n_minmax", std::vector{1, 999});
  const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
  const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;

  // For freq-dependent operators:
  if (h->freqDependantQ && !eachFreqQ)
    h->updateFrequency(const_omega);

  // do RPA:
  std::unique_ptr<ExternalField::TDHF> dV{nullptr};
  if (rpaQ) {
    dV = std::make_unique<ExternalField::TDHF>(h.get(), wf.getHF());
    if (!eachFreqQ)
      dV->solve_core(const_omega);
  }

  std::cout << "\nStructure radiation and normalisation of states:\n";
  if (n_min > 1)
    std::cout << "Including from n = " << n_min << "\n";
  if (n_max < 999)
    std::cout << "Including to n = " << n_max << "\n";
  std::cout << "h=" << h->name() << " (reduced ME)\n";
  std::cout << "Evaluated at ";
  if (eachFreqQ)
    std::cout << "each transition frequency\n";
  else
    std::cout << "constant frequency: w = " << const_omega << "\n";
  if (spline_legs)
    std::cout << "Using splines for diagram legs (external states)\n";
  else
    std::cout << "Using valence states for diagram legs (external states)\n";

  // Lambda to format the output
  const auto printer = [rpaQ](auto str, auto t, auto dv) {
    if (rpaQ)
      printf(" %6s: %12.5e  %12.5e\n", str, t, dv);
    else
      printf(" %6s: %12.5e\n", str, t);
    std::cout << std::flush;
    // nb: the 'flush' is to force a cout flush; particularly when piping
    // output to a file, this wasn't happening soon enough
  };

  if (wf.core.empty() || wf.valence.empty() || wf.basis.empty())
    return;

  // Find core/valence energy: allows distingush core/valence states
  const auto en_core = wf.en_coreval_gap();

  // ----------- ** Actual Calculations ** -----------

  // Construct SR object:
  MBPT::StructureRad sr(wf.basis, en_core, {n_min, n_max});
  std::cout << std::flush;

  struct Output {
    std::string ab{};
    double t0{0.0};
    double t0dv{0.0};
    double sr{0.0};
    double n{0.0};
    double srdv{0.0};
    double ndv{0.0};
  };
  std::vector<Output> out;

  // Loop through all valence states, calc SR+NS
  for (const auto &v : wf.valence) {
    for (const auto &w : wf.valence) {
      if (h->isZero(w.k, v.k))
        continue;

      if (only_diagonal && w != v)
        continue;
      if (!print_both && v > w)
        continue;

      Output t_out;

      // Option to use splines (or valence states) to compute Struc Rad (use
      // splines for legs)
      const auto ws = std::find(cbegin(wf.basis), cend(wf.basis), w);
      const auto vs = std::find(cbegin(wf.basis), cend(wf.basis), v);
      if (spline_legs && (ws == cend(wf.basis) || vs == cend(wf.basis))) {
        std::cout << "Don't have requested spline for: " << w.symbol() << "-"
                  << v.symbol() << "\n";
        continue;
      }
      const auto *vp = spline_legs ? &*vs : &v;
      const auto *wp = spline_legs ? &*ws : &w;

      IO::ChronoTimer timer("time");

      std::cout << "\n" << h->rme_symbol(w, v) << ":\n";
      t_out.ab = h->rme_symbol(w, v);

      const auto ww = eachFreqQ ? std::abs(wp->en() - vp->en()) : const_omega;
      if (eachFreqQ && h->freqDependantQ) {
        h->updateFrequency(ww);
      }
      if (eachFreqQ && rpaQ) {
        if (dV->get_eps() > 1.0e-3)
          dV->clear();
        dV->solve_core(ww);
      }

      // Zeroth-order MEs:
      const auto twvs = h->reducedME(*ws, *vs); // splines here
      const auto twv = h->reducedME(w, v);
      const auto dvs = dV ? twvs + dV->dV(*wp, *vp) : 0.0;
      const auto dv = dV ? twv + dV->dV(w, v) : 0.0;
      printer("t(spl)", twvs, dvs);
      printer("t(val)", twv, dv);

      t_out.t0 = spline_legs ? twvs : twv;
      t_out.t0dv = spline_legs ? dvs : dv;

      // "Top" + "Bottom" SR terms:
      const auto [tb, tb_dv] = sr.srTB(h.get(), *wp, *vp, ww, dV.get());
      printer("SR(TB)", tb, tb_dv);
      // "Centre" SR term:
      const auto [c, c_dv] = sr.srC(h.get(), *wp, *vp, dV.get());
      printer("SR(C)", c, c_dv);

      std::cout << "========\n";
      printer("SR", tb + c, tb_dv + c_dv);

      t_out.sr = tb + c;
      t_out.srdv = tb_dv + c_dv;

      // "Normalisation"
      const auto [n, n_dv] = sr.norm(h.get(), *wp, *vp, dV.get());
      printer("Norm", n, n_dv);

      t_out.n = n;
      t_out.ndv = n_dv;

      printer("Total", tb + c + n, tb_dv + c_dv + n_dv);
      printer("as %", 100.0 * (tb + c + n) / twvs,
              100.0 * (tb_dv + c_dv + n_dv) / dvs);
      out.push_back(t_out);
    }
  }
  std::cout << "\n";
  std::cout << "Structure Radiation + Normalisation of states.\n"
               "Reduced matrix elements (au)\n"
            << "               t0         SR         Norm     ";
  if (rpaQ)
    std::cout << " |  t0+dV      SR+dV      Norm+dV";
  std::cout << "\n";
  for (auto &[ab, t0, t0dv, sr0, n, srdv, ndv] : out) {
    printf("%10s %10.3e %10.3e %10.3e", ab.c_str(), t0, sr0, n);
    if (rpaQ) {
      printf(" | %10.3e %10.3e %10.3e", t0dv, srdv, ndv);
    }
    std::cout << "\n";
  }

  return;
}

//******************************************************************************
void calculateLifetimes(const IO::InputBlock &input, const Wavefunction &wf) {
  std::cout << "\nLifetimes:\n";

  input.checkBlock_old({"E1", "E2", "rpa", "StrucRadNorm"});
  const auto doE1 = input.get("E1", true);
  const auto doE2 = input.get("E2", false);
  const auto rpaQ = input.get("rpa", true);
  if (doE1 && !doE2)
    std::cout << "Including E1 only.\n";
  if (!doE1 && doE2)
    std::cout << "Including E2 only.\n";

  DiracOperator::E1 he1(*(wf.rgrid));
  DiracOperator::Ek he2(*(wf.rgrid), 2);
  const auto alpha = wf.alpha;
  const auto alpha3 = alpha * alpha * alpha;
  const auto alpha2 = alpha * alpha;
  auto dVE1 = ExternalField::TDHF(&he1, wf.getHF());
  auto dVE2 = ExternalField::TDHF(&he2, wf.getHF());

  // Construct SR object:
  std::unique_ptr<MBPT::StructureRad> sr(nullptr);
  const auto srQ = input.get("StrucRadNorm", false);
  if (srQ) {
    std::cout << "Including Structure Radiation + Normalisation\n";
    sr = std::make_unique<MBPT::StructureRad>(wf.basis, wf.en_coreval_gap());
  }

  const auto to_s = PhysConst::time_s;

  struct Data {
    Data(const std::string &s, double t) : state(s), tau(t){};
    std::string state;
    double tau;
  };
  std::vector<Data> data;
  for (const auto &Fa : wf.valence) {
    std::cout << "\n" << Fa.symbol() << "\n";
    auto Gamma = 0.0;

    if (doE1) {
      for (const auto &Fn : wf.valence) {
        if (Fn.en() >= Fa.en() || he1.isZero(Fn.k, Fa.k))
          continue;
        const auto w = Fa.en() - Fn.en();
        if (rpaQ)
          dVE1.solve_core(w, 40);
        auto d = he1.reducedME(Fn, Fa) + dVE1.dV(Fn, Fa);
        if (sr) {
          // include SR.
          const auto [tb, tbx] = sr->srTB(&he1, Fn, Fa);
          const auto [c, cx] = sr->srC(&he1, Fn, Fa);
          const auto [n, nx] = sr->norm(&he1, Fn, Fa);
          d += (tb + c + n);
        }
        const auto g_n = (4.0 / 3) * w * w * w * d * d / (Fa.twojp1());
        Gamma += g_n;
        std::cout << "  E1 --> " << Fn.symbol() << ": ";
        printf("w=%7.5f, |d|=%7.5f, g=%10.4eau\n", w, std::abs(d),
               g_n * alpha3);
      }
    }
    if (doE2) {
      for (const auto &Fn : wf.valence) {
        if (Fn.en() >= Fa.en() || he2.isZero(Fn.k, Fa.k))
          continue;
        const auto w = Fa.en() - Fn.en();
        if (rpaQ)
          dVE2.solve_core(w, 40);
        const auto d = he2.reducedME(Fn, Fa) + dVE2.dV(Fn, Fa);
        const auto g_n = (1.0 / 15) * w * w * w * w * w * d * d / (Fa.twojp1());
        Gamma += g_n * alpha2;
        std::cout << "  E2 --> " << Fn.symbol() << ": ";
        printf("w=%7.5f, |q|=%7.5f, g=%10.4eau\n", w, std::abs(d),
               g_n * alpha3 * alpha2);
      }
    }

    printf("Gamma = %10.4eau = %10.4e/s\n", Gamma * alpha3,
           Gamma * alpha3 / to_s);
    const auto tau = to_s / alpha3 / Gamma;
    printf("tau = %10.4es\n", tau);
    data.emplace_back(Fa.symbol(true), tau * 1.0e9);
  }
  std::cout << "\nLifetimes (summary), in ns:\n";
  for (const auto &[s, t] : data) {
    printf(" %9s  %10.4e\n", s.c_str(), t);
  }
  std::cout << "\n";
}

//******************************************************************************
//******************************************************************************
std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(const IO::InputBlock &input, const Wavefunction &wf,
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
generateOperator(std::string_view oper_name, const IO::InputBlock &input,
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

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_E1(const IO::InputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock_old({"gauge"});
  auto gauge = input.get<std::string>("gauge", "lform");
  if (gauge != "vform")
    return std::make_unique<E1>(*(wf.rgrid));
  // std::cout << "(v-form [velocity gauge])\n";
  return std::make_unique<E1v>(wf.alpha, 0.0);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_Ek(const IO::InputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock_old({"k"});
  auto k = input.get("k", 1);
  return std::make_unique<Ek>(*(wf.rgrid), k);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_M1(const IO::InputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock_old({});
  return std::make_unique<M1>(*(wf.rgrid), wf.alpha, 0.0);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_hfsA(const IO::InputBlock &input, const Wavefunction &wf, bool print) {
  using namespace DiracOperator;
  input.checkBlock_old({"mu", "I", "rrms", "F(r)", "parity", "l", "gl", "mu1",
                        "gl1", "l1", "l2", "I1", "I2", "printF", "screening"});
  auto isotope = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  auto mu = input.get("mu", isotope.mu);
  auto I_nuc = input.get("I", isotope.I_N);
  auto r_rmsfm = input.get("rrms", wf.get_rrms());
  auto r_nucfm = std::sqrt(5.0 / 3) * r_rmsfm;
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
  if (Fr_str == "ball") {
    Fr = Hyperfine::sphericalBall_F();
  } else if (Fr_str == "shell") {
    Fr = Hyperfine::sphericalShell_F();
  } else if (Fr_str == "pointlike" || Fr_str == "point") {
    Fr = Hyperfine::pointlike_F();
  } else if (Fr_str == "VolotkaBW") {
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

    auto mu1 = input.get<double>("mu1", 1.0);
    auto gl1 = input.get<int>("gl1", -1); // 1 or 0 (p or n)
    if (gl1 != 0 && gl1 != 1) {
      std::cout << "FAIL: in " << input.name() << " " << Fr_str
                << "; have gl1=" << gl1 << " but need 1 or 0\n";
      return std::make_unique<NullOperator>(NullOperator());
    }
    auto l1 = input.get<double>("l1", -1.0);
    auto l2 = input.get<double>("l2", -1.0);
    auto I1 = input.get<double>("I1", -1.0);
    auto I2 = input.get<double>("I2", -1.0);

    Fr = Hyperfine::doublyOddBW_F(mu, I_nuc, mu1, I1, l1, gl1, I2, l2);
  } else {
    // read from a file
    const auto [rin, F_of_rin] = IO::FRW::readFile_xy_PoV(Fr_str);
    // interpolate F(r) onto our grid:
    const auto F_r = Interpolator::interpolate(rin, F_of_rin, wf.rgrid->r());
    if (F_r.empty())
      return std::make_unique<NullOperator>();
    return std::make_unique<HyperfineA>(mu, I_nuc, r_nucau, *(wf.rgrid), F_r);
  }

  auto print_FQ = input.get<bool>("printF", false);
  if (print_FQ) {
    std::ofstream of(Fr_str + ".txt");
    for (auto r : wf.rgrid->r()) {
      of << r * PhysConst::aB_fm << " "
         << Fr(r * PhysConst::aB_fm, r_nucau * PhysConst::aB_fm) << "\n";
    }
  }

  return std::make_unique<HyperfineA>(mu, I_nuc, r_nucau, *(wf.rgrid), Fr);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_hfsK(const IO::InputBlock &input, const Wavefunction &wf, bool print) {
  using namespace DiracOperator;

  bool ok = true;
  input.checkBlock_old({"K", "rrms", "F(r)", "gQ"});
  // gQ is g-factor (for magnetic), quadrupole moment for electric...
  const auto k = input.get("K", 0);
  if (k == 0) {
    std::cout
        << "\nFAIL 602: Bad K in hfsK: must be >0 (required: 1=hfsA, 2=hfsB)\n";
    ok = false;
  }

  const auto gQ = input.get("gQ", 1.0);
  if (print) {
    std::cout << "\nHyperfine k=" << k;
    if (k % 2 == 0) {
      std::cout << " (electric) ";
    } else {
      std::cout << " (magnetic) ";
    }
    std::cout << "w/ gQ = " << gQ << "\n";
  }

  const auto isotope = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  const auto r_rmsfm = input.get("rrms", isotope.r_rms);
  const auto r_nucfm = std::sqrt(5.0 / 3) * r_rmsfm;
  const auto r_nucau = r_nucfm / PhysConst::aB_fm;

  auto Fr_str = input.get<std::string>("F(r)", "ball");
  const auto Fr = (Fr_str == "pointlike" || Fr_str == "point") ?
                      Hyperfine::pointlike_F() :
                      Fr_str == "shell" ?
                      Hyperfine::sphericalShell_F() :
                      Fr_str == "ball" ? Hyperfine::sphericalBall_F() :
                                         Hyperfine::pointlike_F();
  if (Fr_str != "ball" && Fr_str != "shell")
    Fr_str = "pointlike";
  else
    std::cout << "Warning: F(r) not correct for k!=1 XXX \n";
  if (print) {
    std::cout << "w/ " << Fr_str << " for F(r), r_N = " << r_nucfm << "fm "
              << " (rrms=" << r_rmsfm << "fm)\n";
  }

  if (!ok)
    return std::make_unique<NullOperator>();

  return std::make_unique<HyperfineK>(k, gQ, r_nucau, *wf.rgrid, Fr);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_r(const IO::InputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock_old({"power"});
  auto power = input.get("power", 1.0);
  std::cout << "r^(" << power << ")\n";
  return std::make_unique<RadialF>(*(wf.rgrid), power);
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_pnc(const IO::InputBlock &input, const Wavefunction &wf, bool) {
  using namespace DiracOperator;
  input.checkBlock_old({"c", "t"});
  const auto r_rms = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
  const auto c = input.get("c", Nuclear::c_hdr_formula_rrms_t(r_rms));
  const auto t = input.get("t", Nuclear::default_t);
  return std::make_unique<PNCnsi>(c, t, *(wf.rgrid), 1.0, "iQwe-11");
}

//------------------------------------------------------------------------------
std::unique_ptr<DiracOperator::TensorOperator>
generate_Hrad(const IO::InputBlock & /*input*/, const Wavefunction & /*wf*/,
              bool) {
  std::cout << "\nFAIL:: generate_Hrad() need implementing!\n";
  return nullptr;
  // using namespace DiracOperator;
  // input.checkBlock_old(
  //     {"Simple", "Ueh", "SE_h", "SE_l", "SE_m", "rcut", "scale_rN"});
  // const auto x_Simple = input.get("Simple", 0.0);
  // const auto x_Ueh = input.get("Ueh", 1.0);
  // const auto x_SEe_h = input.get("SE_h", 1.0);
  // const auto x_SEe_l = input.get("SE_l", 1.0);
  // const auto x_SEm = input.get("SE_m", 1.0);
  // const auto rcut = input.get("rcut", 5.0);
  // const auto scale_rN = input.get("scale_rN", 1.0);
  // const auto r_rms_Fermi = scale_rN * wf.get_nuclearParameters().r_rms;
  // const auto Hel = RadiativePotential::form_Hel(wf.rgrid->r() , x_Simple,
  // x_Ueh,
  //                                               x_SEe_h, x_SEe_l,
  //                                               r_rms_Fermi, wf.Znuc(),
  //                                               wf.alpha, rcut);
  // const auto Hmag = RadiativePotential::form_Hmag(
  //     wf.rgrid->r() , x_SEm, r_rms_Fermi, wf.Znuc(), wf.alpha, rcut);
  // return std::make_unique<Hrad>(Hel, Hmag);
  // // return std::make_unique<Hrad_el>(Hel);
}

//******************************************************************************

} // namespace Module
