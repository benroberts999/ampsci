#include "Modules/Module_matrixElements.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "DiracOperator/Operators.hpp"
#include "HF/ExternalField.hpp"
#include "IO/UserInput.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>

namespace Module {

//******************************************************************************
void matrixElements(const UserInputBlock &input, const Wavefunction &wf) {

  std::string ThisModule = "MatrixElements::";
  auto operator_str = input.name().substr(ThisModule.length());
  // nb: "check" is done in 'generate operator'

  if (operator_str == "SecondOrder" || operator_str == "Sigma2") {
    // Special case
    return SecondOrder(input, wf);
  }

  const bool radial_int = input.get("radialIntegral", false);

  // spacial case: HFS A (MHz)
  bool AhfsQ = (operator_str == "hfs" && !radial_int);

  auto which_str =
      radial_int ? "(radial integral). " : AhfsQ ? " A (MHz). " : "(reduced). ";

  const auto h = generateOperator(operator_str, input, wf);
  std::cout << "\n"
            << ThisModule << which_str << " Operator: " << h->name() << "\n";
  std::cout << "Units: " << h->units() << "\n";

  const bool print_both = input.get("printBoth", false);
  const bool diagonal_only = input.get("onlyDiagonal", false);

  const bool rpaQ = input.get("rpa", true);
  const auto omega = input.get("omega", 0.0);
  const bool eachFreqQ = omega < 0.0;

  const auto factor = input.get("factor", 1.0);
  if (factor != 1.0)
    std::cout << "With extra factor: " << factor << "\n";

  // const auto units = input.get<std::string>("units", "au"); // can remove?
  // const auto unit = (units == "MHz" || AhfsQ) ? PhysConst::Hartree_MHz : 1.0;

  auto rpa =
      HF::ExternalField(h.get(), wf.core_orbitals,
                        NumCalc::add_vectors(wf.vnuc, wf.vdir), wf.get_alpha());
  std::unique_ptr<HF::ExternalField> rpa0; // for first-order

  if (!eachFreqQ && rpaQ) {
    rpa.solve_TDHFcore(omega, 1, false);
    rpa0 =
        std::make_unique<HF::ExternalField>(rpa); // store first-order snapshot
    rpa.solve_TDHFcore(omega);
  } else {
    rpa0 = std::make_unique<HF::ExternalField>(rpa); // Solved later
  }

  // Fb -> Fa = <a||h||b>
  for (const auto &Fb : wf.valence_orbitals) {
    for (const auto &Fa : wf.valence_orbitals) {
      if (h->isZero(Fa.k, Fb.k))
        continue;
      if (diagonal_only && Fb != Fa)
        continue;
      if (!print_both && Fb > Fa)
        continue;
      if (eachFreqQ && rpaQ) {
        auto w = std::abs(Fa.en - Fb.en);
        rpa0->clear_dPsi();
        rpa0->solve_TDHFcore(w, 1, false); // wastes a little time
        rpa.solve_TDHFcore(w); // re-solve at new frequency (not from scratch)
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
void calculateBohrWeisskopf(const UserInputBlock &input, const Wavefunction &wf)
//
{
  using namespace DiracOperator;

  UserInputBlock point_in("MatrixElements::hfs", input);
  UserInputBlock ball_in("MatrixElements::hfs", input);
  UserInputBlock BW_in("MatrixElements::hfs", input);
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
  for (const auto &phi : wf.valence_orbitals) {
    auto Ap = Hyperfine::hfsA(hp.get(), phi);
    auto Ab = Hyperfine::hfsA(hb.get(), phi);
    auto Aw = Hyperfine::hfsA(hw.get(), phi);
    auto Fball = ((Ab / Ap) - 1.0) * 100.0; //* M_PI * PhysConst::c;
    auto Fbw = ((Aw / Ap) - 1.0) * 100.0;   //* M_PI * PhysConst::c;
    // printf("%6s: %9.1f  %9.1f  %9.1f | %8.4f  %8.4f   %9.6f\n",
    printf("%7s| %12.5e %12.5e %12.5e | %9.5f  %9.5f \n", phi.symbol().c_str(),
           Ap, Ab, Aw, Fball, Fbw);
  }
}

//******************************************************************************
void SecondOrder(const UserInputBlock &input, const Wavefunction &wf) {

  input.checkBlock({"lmax", "kmax"});
  auto kmax = input.get("kmax", 10);
  auto lmax = input.get("lmax", 10);

  // std::vector<DiracSpinor> core;
  // for (const auto &Fb : wf.basis) {
  //   if (wf.isInCore(Fb))
  //     core.push_back(Fb);
  // }
  const auto &core = wf.core_orbitals;

  std::cout << "\nMBPT(2): Valence energy shifts.\n";
  std::cout << "Matrix elements <v|Sigma(2)|v>:\n";

  if (wf.basis.empty())
    std::cout << "FAIL 125 in Module::SecondOrder: There is no basis! - I need "
                 "a basis to calculate MBPT. Try again.\n";

  // auto Yk = Coulomb::YkTable(&wf.rgrid, &wf.basis, &core);

  for (const auto &v : wf.valence_orbitals) {
    if (v.l() > lmax)
      continue;

    std::vector<double> delta_b(core.size());
#pragma omp parallel for
    for (auto ib = 0ul; ib < core.size(); ib++) {
      double sigma_b = 0.0;
      for (int k = 0; k <= kmax; ++k) {
        double sigma_k = 0.0;
        auto f = (2 * k + 1) * v.twojp1();
        const auto &b = core[ib];
        for (const auto &n : wf.basis) {
          if (wf.isInCore(n))
            continue;
          for (const auto &a : core) {
            const auto zx = Coulomb::Wk_abcd(v, n, a, b, k) *
                            Coulomb::Qk_abcd(v, n, a, b, k);
            // const auto [kmin_nb, kmax_nb] = Yk.k_minmax(n, b);
            // if (k < kmin_nb || k > kmax_nb)
            //   continue;
            // const auto &yknb = Yk.get_yk_ab(k, n, b);
            // const auto &yna = Yk.get_y_ab(n, a);
            // const auto zx = Coulomb::Wk_abcd(v, n, a, b, k, yknb, yna) *
            //                 Coulomb::Qk_abcd(v, n, a, b, k, yknb);
            const auto dele = v.en + n.en - a.en - b.en;
            sigma_k += zx / dele;
          } // a
          for (const auto &m : wf.basis) {
            if (wf.isInCore(m))
              continue;
            const auto zx = Coulomb::Wk_abcd(v, b, m, n, k) *
                            Coulomb::Qk_abcd(v, b, m, n, k);
            // const auto [kmin_nb, kmax_nb] = Yk.k_minmax(n, b);
            // if (k < kmin_nb || k > kmax_nb)
            //   continue;
            // const auto &ykbn = Yk.get_yk_ab(k, n, b);
            // const auto &ybm = Yk.get_y_ab(m, b);
            // const auto zx = Coulomb::Wk_abcd(v, b, m, n, k, ykbn, ybm) *
            //                 Coulomb::Qk_abcd(v, b, m, n, k, ykbn);
            const auto dele = m.en + n.en - v.en - b.en;
            sigma_k -= zx / dele;
          } // m
        }   // n
        sigma_b += sigma_k / f;
      } // k
      delta_b[ib] = sigma_b;
    } // b

    auto delta = std::accumulate(delta_b.begin(), delta_b.end(), 0.0);

    printf("%7s| %9.6f %+9.6f = %9.6f = %9.2f\n", v.symbol().c_str(), v.en,
           delta, (v.en + delta), (v.en + delta) * PhysConst::Hartree_invcm);
    if (std::abs(delta / v.en) > 0.5)
      std::cout << "      *** Warning: delta V. large?\n";
  }
}

//******************************************************************************
std::unique_ptr<DiracOperator::TensorOperator>
generateOperator(const std::string &operator_str, const UserInputBlock &input,
                 const Wavefunction &wf) {
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
      return std::make_unique<E1>(E1(wf.rgrid));
    // std::cout << "(v-form [velocity gauge])\n";
    return std::make_unique<E1_vform>(E1_vform(wf.rgrid, wf.get_alpha()));
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
    std::cout << "Sorry, check back soon for M1 :(\n";
    // return std::make_unique<M1Operator>(M1Operator());
    return std::make_unique<NullOperator>(NullOperator());
  }

  std::cerr << "\nFAILED to find operator: " << ThisModule
            << " in generateOperator. Returning NULL operator (0)\n";
  return std::make_unique<NullOperator>(NullOperator());
}

//******************************************************************************
void calculateLifetimes(const UserInputBlock &input, const Wavefunction &wf) {
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
  auto alpha = wf.get_alpha();
  auto alpha3 = alpha * alpha * alpha;
  auto alpha2 = alpha * alpha;
  auto dVE1 = HF::ExternalField(&he1, wf.core_orbitals,
                                NumCalc::add_vectors(wf.vnuc, wf.vdir), alpha);
  auto dVE2 = HF::ExternalField(&he2, wf.core_orbitals,
                                NumCalc::add_vectors(wf.vnuc, wf.vdir), alpha);

  auto to_s = PhysConst::time_s;

  for (const auto &Fa : wf.valence_orbitals) {
    std::cout << "\n" << Fa.symbol() << "\n";
    auto Gamma = 0.0;

    if (doE1) {
      for (const auto &Fn : wf.valence_orbitals) {
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
      for (const auto &Fn : wf.valence_orbitals) {
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
