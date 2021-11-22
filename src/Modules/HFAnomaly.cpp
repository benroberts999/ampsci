#include "Modules/HFAnomaly.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/matrixElements.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace Module {

static void calc_thing(const DiracSpinor &Fv, double e_targ, double r0,
                       double mu, double I, int l, int gl);

//******************************************************************************
void HFAnomaly(const IO::InputBlock &input, const Wavefunction &wf) {

  input.checkBlock_old({"rpa", "options", "A"});

  const auto rpa = input.get("rpa", false);
  const auto Alist = input.get("A", std::vector<int>{});

  auto options = input.getBlock("options");
  auto sub_input = IO::InputBlock("hfs", {});
  if (options) {
    sub_input.add(options->options());
  }

  auto point_in = sub_input;
  auto ball_in = sub_input;
  auto sp_in = sub_input;

  point_in.add("F(r)=pointlike;");
  ball_in.add("F(r)=ball;");
  const auto doubly_odd = (wf.Anuc() % 2 == 0);
  if (doubly_odd)
    sp_in.add("F(r)=doublyOddBW;");
  else
    sp_in.add("F(r)=VolotkaBW;");

  // Uses generateOperator (from Module::MatrixElements)
  // --> Probably just as easy to do from scratch?
  const auto hpt = generateOperator(point_in, wf, false);
  const auto hbl = generateOperator(ball_in, wf, false);
  const auto hsp = generateOperator(sp_in, wf, true);

  // Struct to store hfs constant for the reference isotope (point, ball,
  // single-particle)
  struct ThreeAs {
    ThreeAs(double a, double b, double c) : pt(a), bl(b), sp(c){};
    double pt, bl, sp;
  };

  // calc mu, I, gI
  const auto isotope0 = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  const auto mu0 = sub_input.get("mu", isotope0.mu);
  const auto I0 = sub_input.get("I", isotope0.I_N);
  const auto r_rms0 = sub_input.get("rrms", isotope0.r_rms);
  const auto gI0 = mu0 / I0;

  // RPA:
  std::unique_ptr<ExternalField::DiagramRPA> rpap{nullptr}, rpab{nullptr},
      rpas{nullptr};
  if (rpa) {
    std::cout << "\nIncluding RPA (diagram method) - must have basis\n";
    rpap = std::make_unique<ExternalField::DiagramRPA>(hpt.get(), wf.basis,
                                                       wf.core, wf.identity());
    rpab = std::make_unique<ExternalField::DiagramRPA>(hbl.get(), rpap.get());
    rpas = std::make_unique<ExternalField::DiagramRPA>(hsp.get(), rpap.get());
    std::cout << "Solving RPA core for point, ball, SP:\n";
    rpap->solve_core(0.0, 100, true);
    rpab->solve_core(0.0, 100, true);
    rpas->solve_core(0.0, 100, true);
  }

  // Calculate HFS for each state of reference isotope: including BW effect
  std::vector<ThreeAs> As;
  std::cout << "\nReference: A = " << wf.Anuc() << ", gI = " << gI0
            << ", r_rms = " << r_rms0 << "\n\n";
  std::cout
      << "      A(point)     A(ball)      A(SP)       | e(ball)  e(SP) [%]\n";
  for (const auto &Fv : wf.valence) {
    auto point = DiracOperator::HyperfineA::hfsA(hpt.get(), Fv);
    auto ball = DiracOperator::HyperfineA::hfsA(hbl.get(), Fv);
    auto sp = DiracOperator::HyperfineA::hfsA(hsp.get(), Fv);
    if (rpa) {
      const auto a = DiracOperator::HyperfineA::convertRMEtoA(Fv, Fv);
      point += a * rpap->dV(Fv, Fv);
      ball += a * rpab->dV(Fv, Fv);
      sp += a * rpas->dV(Fv, Fv);
    }
    As.emplace_back(point, ball, sp);
    const auto eps_ball = 100.0 * (ball - point) / point;
    const auto eps_sp = 100.0 * (sp - point) / point;
    printf("%-4s %12.5e %12.5e %12.5e | %7.4f  %7.4f\n",
           Fv.shortSymbol().c_str(), point, ball, sp, eps_ball, eps_sp);
  }

  // For each known odd isotope, calculate HFS and HFS-anomaly
  // Note: Only odd, since then can simply guess VolotkaBW
  // For even: choose as the reference isotope
  const auto isotopes = Nuclear::findIsotopeList(wf.Znuc());
  for (const auto &nuc : isotopes) {
    // only do odd isotopes (even must be reference!)
    if (nuc.A % 2 == 0)
      continue;
    if (!Alist.empty()) {
      // if Alist is given, only run for those isotopes
      if (std::find(Alist.cbegin(), Alist.cend(), nuc.A) == Alist.cend())
        continue;
    }
    auto wfA = Wavefunction(wf.rgrid->params(), {nuc.Z, nuc.A},
                            wf.alpha / PhysConst::alpha);
    const auto gI = nuc.mu / nuc.I_N;
    std::cout << "\nA = " << nuc.A << ", gI = " << gI
              << ", r_rms = " << nuc.r_rms << "\n";
    wfA.solve_core("HartreeFock", 0.0, wf.coreConfiguration_nice());
    wfA.solve_valence(DiracSpinor::state_config(wf.valence));
    wfA.basis = wf.basis; // OK??
    //  wfA.formBasis({"50spd30f", 60, 7, 0.0, 0.0, 30.0, false});

    const auto hpt2 = generateOperator({"hfs", "F(r)=pointlike;"}, wfA, false);
    const auto hbl2 = generateOperator({"hfs", "F(r)=ball;"}, wfA, false);
    const auto hsp2 = generateOperator({"hfs", "F(r)=VolotkaBW;"}, wfA, false);

    std::unique_ptr<ExternalField::DiagramRPA> rpap2{nullptr}, rpab2{nullptr},
        rpas2{nullptr};
    if (rpa) {
      std::cout << "Including RPA (diagram method) - must have basis\n";
      // OK ? or need run with new basis!?
      rpap2 = std::make_unique<ExternalField::DiagramRPA>(
          hpt2.get(), wfA.basis, wfA.core, wfA.identity());
      rpab2 =
          std::make_unique<ExternalField::DiagramRPA>(hbl2.get(), rpap.get());
      rpas2 =
          std::make_unique<ExternalField::DiagramRPA>(hsp2.get(), rpap.get());
      std::cout << "Solving RPA core for point, ball, SP:\n";
      rpap2->solve_core(0.0, 100, true);
      rpab2->solve_core(0.0, 100, true);
      rpas2->solve_core(0.0, 100, true);
    }

    // for (const auto &Fv : wf.valence)
    std::cout << "\nA = " << nuc.A << ", gI = " << gI
              << ", r_rms = " << nuc.r_rms << "\n";
    std::cout << "      A(point)      e(ball)   e(SP)   | 1D2(ball) "
                 "1D2(SP) [%]\n";
    for (std::size_t i = 0; i < wf.valence.size(); ++i) {
      // Make sure the valence states are in correct order!
      const auto &Fv = *wfA.getState(wf.valence[i].shortSymbol());
      const auto [pt0, bl0, sp0] = As[i];
      auto point = DiracOperator::HyperfineA::hfsA(hpt2.get(), Fv);
      auto ball = DiracOperator::HyperfineA::hfsA(hbl2.get(), Fv);
      auto sp = DiracOperator::HyperfineA::hfsA(hsp2.get(), Fv);
      if (rpa) {
        const auto a = DiracOperator::HyperfineA::convertRMEtoA(Fv, Fv);
        point += a * rpap2->dV(Fv, Fv);
        ball += a * rpab2->dV(Fv, Fv);
        sp += a * rpas2->dV(Fv, Fv);
      }
      // As.emplace_back(point, ball, sp);
      const auto eps_ball = 100.0 * (ball - point) / point;
      const auto eps_sp = 100.0 * (sp - point) / point;
      const auto D12_bl = 100.0 * ((bl0 / ball) * (gI / gI0) - 1.0);
      const auto D12_sp = 100.0 * ((sp0 / sp) * (gI / gI0) - 1.0);
      if (std::abs(eps_ball) > 0.001)
        printf("%-4s %12.5e | %8.5f %8.5f | %8.5f  %8.5f\n",
               Fv.shortSymbol().c_str(), point, eps_ball, eps_sp, D12_bl,
               D12_sp);
      else
        printf("%-4s %12.5e | %8.1e %8.1e | %8.1e  %8.1e\n",
               Fv.shortSymbol().c_str(), point, eps_ball, eps_sp, D12_bl,
               D12_sp);
    }
  }
}

//******************************************************************************
void HF_rmag(const IO::InputBlock &input, const Wavefunction &wf) {
  // For isotope 1 and 2
  // Loops over many values for magnetic radius of isotope 1, Rmag(1).
  // For each, finds Rmag(2) that reproduces a given hyperfine anomaly.
  // Potential to use two states to find rmag. If not, still interesting?

  std::cout << "\nTuning Rmag to fit hyperfine anomaly\n";

  input.checkBlock_old({"n", "kappa", "A2", "1D2", "rpa", "num_steps", "mu1",
                        "mu2", "I1", "I2", "eps_targ", "e1", "e2"});

  // A(1) is wf
  // A(2) is wf2
  const auto A2 = input.get<int>("A2", wf.Anuc());
  const auto eps_t = input.get("eps_targ", 1.0e-4);

  const auto n = input.get("n", 1);
  const auto kappa = input.get("kappa", -1);

  auto wf2 = Wavefunction(
      wf.rgrid->params(),
      {wf.Znuc(), A2, Nuclear::parseType(wf.get_nuclearParameters().type)},
      wf.alpha / PhysConst::alpha);
  wf2.solve_core("HartreeFock", 0.0, wf.coreConfiguration_nice());
  wf2.solve_valence(DiracSpinor::state_config(wf.valence));
  wf2.basis = wf.basis; // OK??
  std::cout << "A1 = " << wf.nuclearParams() << "\n";
  std::cout << "A2 = " << wf2.nuclearParams() << "\n";

  const auto iso_1 = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  const auto iso_2 = Nuclear::findIsotopeData(wf2.Znuc(), wf2.Anuc());

  // calc mu, I, gI
  const auto mu1 = input.get("mu1", iso_1.mu);
  const auto mu2 = input.get("mu2", iso_2.mu);
  const auto I1 = input.get("I1", iso_1.I_N);
  const auto I2 = input.get("I2", iso_2.I_N);
  std::cout << "mu1=" << mu1 << ", mu1=" << mu2 << "\n";

  const auto g1 = mu1 / I1;
  const auto g2 = mu2 / I2;

  const auto r0_fm = input.get("rrms", iso_1.r_rms);
  const auto convert = std::sqrt(5.0 / 3) / PhysConst::aB_fm;
  const auto rN0_au = r0_fm * convert;
  const auto rN2 = iso_2.r_rms * convert;

  // Parameters for Single-Particle BW effect:
  const auto gl = wf.Znuc() % 2 == 0 ? 0 : 1; // unparied proton?
  const auto pi1 = iso_1.parity;
  const auto l1_tmp = int(I1 + 0.5 + 0.0001);
  const auto l1 = ((l1_tmp % 2 == 0) == (pi1 == 1)) ? l1_tmp : l1_tmp - 1;
  const auto pi2 = iso_2.parity;
  const auto l2_tmp = int(I2 + 0.5 + 0.0001);
  const auto l2 = ((l2_tmp % 2 == 0) == (pi2 == 1)) ? l2_tmp : l2_tmp - 1;

  const auto Fr1 = DiracOperator::Hyperfine::volotkaBW_F(mu1, I1, l1, gl);
  const auto Fr2 = DiracOperator::Hyperfine::volotkaBW_F(mu2, I2, l2, gl);
  // const auto Fr1 = DiracOperator::Hyperfine::sphericalBall_F();
  // const auto Fr2 = DiracOperator::Hyperfine::sphericalBall_F();

  const auto &F1v = *wf.getState(n, kappa);
  const auto &F2v = *wf2.getState(n, kappa);

  // nb: can only do diagram RPA for hfs
  const auto rpa = input.get("rpa", false) || input.get("rpa_diagram", false);

  // nb: this is very ineficient...
  std::unique_ptr<ExternalField::DiagramRPA> rpa01{nullptr}, rpa02{nullptr},
      rpa1{nullptr}, rpa2{nullptr}, rpa2a{nullptr}, rpa2b{nullptr};

  // const double x = 0.05;
  // const auto num_steps = input.get("num_steps", 40) + 1;
  // const auto dr = 5.0 * x * rN0_au / (num_steps + 1);

  const auto h01 = DiracOperator::HyperfineA(
      mu1, I1, 0.0, *wf.rgrid, DiracOperator::Hyperfine::pointlike_F());
  const auto h02 = DiracOperator::HyperfineA(
      mu2, I2, 0.0, *wf.rgrid, DiracOperator::Hyperfine::pointlike_F());
  double dv01 = 0.0;
  double dv02 = 0.0;
  if (rpa) {
    rpa01 = std::make_unique<ExternalField::DiagramRPA>(&h01, wf.basis, wf.core,
                                                        wf.identity());
    rpa02 = std::make_unique<ExternalField::DiagramRPA>(&h02, wf.basis, wf.core,
                                                        wf.identity());
    rpa01->solve_core(0.0);
    rpa02->solve_core(0.0);
    auto a = DiracOperator::HyperfineA::convertRMEtoA(F1v, F1v);
    dv01 = a * rpa01->dV(F1v, F1v);
    dv02 = a * rpa02->dV(F1v, F1v);
  }
  double a0_1 = h01.hfsA(F1v) + dv01;
  double a0_2 = h02.hfsA(F1v) + dv02;

  // auto rat_min = 1.05;
  // auto rat_max = 1.12;
  auto rat_min = 1.065;
  auto rat_max = 1.095;
  const auto num_steps = input.get("num_steps", 40);
  const auto dr = (rat_max - rat_min) / (num_steps - 1);

  std::cout << "Rc = Sqrt[5/3]r_rms:\n";
  std::cout << "Rc(1) = " << rN0_au * PhysConst::aB_fm << "\n";
  std::cout << "Rc(2) = " << rN2 * PhysConst::aB_fm << "\n";

  std::cout << "\nRunning for " << F1v.symbol() << "\n";
  const auto d12_targ = input.get<double>("1D2", 0.0);
  if (d12_targ != 0.0) {
    std::cout << "A0(1) = " << a0_1 << "\n";
    std::cout << "A0(2) = " << a0_2 << "\n";
    std::cout << "1D2 target: " << d12_targ << "\n";
    std::cout << "Rm/Rc(1) Rm/Rc(2) e1     e2     1D2    eps(D)   del(Rm/Rc)\n";
    bool done_one_yet = false;
    std::cout << std::flush;

    for (double rat = rat_min; rat < rat_max + 0.5 * dr; rat += dr) {
      double rN = rat * rN0_au;
      std::cout << std::flush;

      if (!done_one_yet && rN >= rN0_au) {
        // ensure we do exactly Rm/Rc = 1
        rN = rN0_au;
        done_one_yet = true;
        rat -= dr;
      }

      const auto h1 = DiracOperator::HyperfineA(mu1, I1, rN, *wf.rgrid, Fr1);
      double dv1 = 0.0;
      if (rpa) {
        rpa1 = std::make_unique<ExternalField::DiagramRPA>(&h1, rpa01.get());
        rpa1->grab_tam(rpa01.get()); // don't start from scratch
        rpa1->solve_core(0.0, 100, false);
        auto a = DiracOperator::HyperfineA::convertRMEtoA(F1v, F1v);
        dv1 = a * rpa1->dV(F1v, F1v);
      }

      const auto a1 = h1.hfsA(F1v) + dv1;

      // Use halving-interval method to find rmag(2) that reproduces 1D2,
      // given rmag(1)
      auto r2a = rat_min * rN2;
      auto r2 = rN2;
      auto r2b = rat_max * rN2;
      int tries = 0;
      double eps;
      double del_r;
      double d12;
      double bw1, bw2;
      while (tries++ < 40) {

        // Note: For RPA, could make this 3x faster!
        // Only need to change the 'outer' guess each time, do not need to
        // re-evaluate RPA for the other two guesses! Complicated though??
        const auto h2a =
            DiracOperator::HyperfineA(mu2, I2, r2a, *wf.rgrid, Fr2);
        const auto h2 = DiracOperator::HyperfineA(mu2, I2, r2, *wf.rgrid, Fr2);
        const auto h2b =
            DiracOperator::HyperfineA(mu2, I2, r2b, *wf.rgrid, Fr2);

        double dv2 = 0.0;
        double dv2a = 0.0;
        double dv2b = 0.0;
        if (rpa) {
          // this is very ineficient.. should only need to solve once per run..?
          rpa2a = std::make_unique<ExternalField::DiagramRPA>(&h2a, rpa1.get());
          rpa2 = std::make_unique<ExternalField::DiagramRPA>(&h2, rpa1.get());
          rpa2b = std::make_unique<ExternalField::DiagramRPA>(&h2b, rpa1.get());
          rpa2->grab_tam(rpa02.get());  // don't start from scratch
          rpa2a->grab_tam(rpa02.get()); // don't start from scratch
          rpa2b->grab_tam(rpa02.get()); // don't start from scratch
          rpa2a->solve_core(0.0, 100, false);
          rpa2->solve_core(0.0, 100, false);
          rpa2b->solve_core(0.0, 100, false);
          const auto a = DiracOperator::HyperfineA::convertRMEtoA(F2v, F2v);
          dv2a = a * rpa2a->dV(F2v, F2v);
          dv2 = a * rpa2->dV(F2v, F2v);
          dv2b = a * rpa2b->dV(F2v, F2v);
        }

        const auto a2a = h2a.hfsA(F2v) + dv2a;
        const auto a2 = h2.hfsA(F2v) + dv2;
        const auto a2b = h2b.hfsA(F2v) + dv2b;
        const auto d12a = 100.0 * ((a1 / a2a) * (g2 / g1) - 1.0);
        d12 = 100.0 * ((a1 / a2) * (g2 / g1) - 1.0);
        const auto d12b = 100.0 * ((a1 / a2b) * (g2 / g1) - 1.0);

        bw1 = 100.0 * (a1 - a0_1) / a0_1;
        bw2 = 100.0 * (a2 - a0_2) / a0_2;

        // Halfing interval:
        if ((d12a < d12_targ && d12_targ < d12) ||
            (d12 < d12_targ && d12_targ < d12a)) {
          // r2a = r2a;
          r2b = r2;
          r2 = 0.5 * (r2a + r2b);
        } else if ((d12b < d12_targ && d12_targ < d12) ||
                   (d12 < d12_targ && d12_targ < d12b)) {
          // r2b = r2b;
          r2a = r2;
          r2 = 0.5 * (r2a + r2b);
        } else {
          r2a *= 0.95;
          r2b *= 1.05;
        }

        eps = std::abs((d12 - d12_targ) / d12_targ);
        del_r = 0.5 * std::abs(r2a - r2b);
        if (eps < eps_t) {
          break;
        }
      } // tries
      printf("%.5f  %.5f %6.3f %6.3f %7.4f %.1e  %.1e  [%i]\n", rN / rN0_au,
             r2 / rN2, bw1, bw2, d12, eps, del_r / rN0_au, tries);
    } // rMag
  }

  auto etarg1 = input.get("e1", 0.0);
  if (etarg1 != 0.0) {
    std::cout << "\nFitting rmag to find epsilon for A1: e_targ=" << etarg1
              << "\n";
    calc_thing(F1v, etarg1, rN0_au, mu1, I1, l1, gl);
  }

  auto etarg2 = input.get("e2", 0.0);
  if (etarg2 != 0.0) {
    std::cout << "\nFitting rmag to find epsilon for A2: e_targ=" << etarg2
              << "\n";
    calc_thing(F2v, etarg2, rN2, mu2, I2, l2, gl);
  }
}

//******************************************************************************
static void calc_thing(const DiracSpinor &Fv, double e_targ, double r0,
                       double mu, double I, int l, int gl) {
  // I'm so glad past Ben gave this function such an informative name....
  // I think this finds the r_mag required to reproduce epsilon

  const auto Fr = DiracOperator::Hyperfine::volotkaBW_F(mu, I, l, gl);

  const auto h0 = DiracOperator::HyperfineA(
      mu, I, 0.0, *Fv.rgrid, DiracOperator::Hyperfine::pointlike_F());
  double dv0 = 0.0; // XXX RPA!
  const auto A0 = h0.hfsA(Fv) + dv0;

  std::cout << "nb: no RPA\n";

  int tries = 0.0;
  double r = r0;
  double ra = 0.8 * r0, rb = 1.2 * r0;
  while (tries++ < 100) {
    const auto ha = DiracOperator::HyperfineA(mu, I, ra, *Fv.rgrid, Fr);
    const auto h = DiracOperator::HyperfineA(mu, I, r, *Fv.rgrid, Fr);
    const auto hb = DiracOperator::HyperfineA(mu, I, rb, *Fv.rgrid, Fr);

    double dv = 0.0, dva = 0.0, dvb = 0.0; // XXX Inlcude RPA!
    const auto Aa = ha.hfsA(Fv) + dva;
    const auto A = h.hfsA(Fv) + dv;
    const auto Ab = hb.hfsA(Fv) + dvb;

    auto bwa = 100.0 * (Aa - A0) / A0;
    auto bw = 100.0 * (A - A0) / A0;
    auto bwb = 100.0 * (Ab - A0) / A0;

    // Halfing interval:
    if ((bwa < e_targ && e_targ < bw) || (bw < e_targ && e_targ < bwa)) {
      // ra = ra;
      rb = r;
      r = 0.5 * (ra + rb);
    } else if ((bwb < e_targ && e_targ < bw) || (bw < e_targ && e_targ < bwb)) {
      // rb = rb;
      ra = r;
      r = 0.5 * (ra + rb);
    } else {
      ra *= 0.95;
      rb *= 1.05;
    }

    if (std::abs((bwa - bwb) / bw) < 1.0e-6) {
      // auto c = PhysConst::aB_fm;
      printf(
          "Rm/Rc=%8.6f, e=%6.3f,   del_Rm/Rc=%.1e   eps_r=%.0e  eps_e=%.0e\n",
          r / r0, bw, std::abs(rb - ra) / r0, (rb - ra) / r, (bwa - bwb) / bw);
      break;
    }
  }
}

//******************************************************************************
void calculateBohrWeisskopf(const IO::InputBlock &input,
                            const Wavefunction &wf) {
  using namespace DiracOperator;

  input.checkBlock_old({"rpa", "rpa_diagram", "screening", "hfs_options"});

  // const auto h_options =
  //     IO::InputBlock("hfs_options", input.get<std::string>("hfs_options",
  //     ""));
  // IO::InputBlock point_in("hfs", h_options);
  // IO::InputBlock ball_in("hfs", h_options);
  // IO::InputBlock BW_in("hfs", h_options);
  // point_in.add("F(r)=pointlike");
  // ball_in.add("F(r)=ball");
  // if (wf.Anuc() % 2 == 0)
  //   BW_in.add("F(r)=doublyOddBW");
  // else
  //   BW_in.add("F(r)=VolotkaBW");

  auto options = input.getBlock("hfs_options");
  auto sub_input = IO::InputBlock("hfs", {});
  if (options) {
    sub_input.add(options->options());
  }

  auto point_in = sub_input;
  auto ball_in = sub_input;
  auto BW_in = sub_input;

  point_in.add("F(r)=pointlike;");
  ball_in.add("F(r)=ball;");
  const auto doubly_odd = (wf.Anuc() % 2 == 0);
  if (doubly_odd)
    BW_in.add("F(r)=doublyOddBW;");
  else
    BW_in.add("F(r)=VolotkaBW;");

  auto hp = generateOperator(point_in, wf, false);
  auto hb = generateOperator(ball_in, wf, false);
  auto hw = generateOperator(BW_in, wf);

  // nb: can only do diagram RPA for hfs
  const auto rpa = input.get("rpa", false) || input.get("rpa_diagram", false);

  std::unique_ptr<ExternalField::DiagramRPA> rpap{nullptr}, rpab{nullptr},
      rpaw{nullptr};
  if (rpa) {
    std::cout << "\nIncluding RPA (diagram method) - must have basis\n";
    rpap = std::make_unique<ExternalField::DiagramRPA>(hp.get(), wf.basis,
                                                       wf.core, wf.identity());
    rpab = std::make_unique<ExternalField::DiagramRPA>(hb.get(), rpap.get());
    rpaw = std::make_unique<ExternalField::DiagramRPA>(hw.get(), rpap.get());
    std::cout << "Solving RPA core for point, ball, SP:\n";
    rpap->solve_core(0.0);
    rpab->solve_core(0.0);
    rpaw->solve_core(0.0);
  }

  // only used if calculating screening factors:
  struct ScreenBW {
    ScreenBW(int a, int b1, double c, double d) : n(a), k(b1), b(c), sp(d){};
    int n;
    int k;
    double b;
    double sp;
  };
  std::vector<ScreenBW> bw;

  std::cout << "\nTabulate A (Mhz), and Bohr-Weisskopf effect eps(%): "
            << wf.atom()
            << "\n       |A:      Point         Ball           SP |e:    "
               "Ball         SP\n";
  for (const auto &phi : wf.valence) {
    auto Ap = HyperfineA::hfsA(hp.get(), phi);
    auto Ab = HyperfineA::hfsA(hb.get(), phi);
    auto Aw = HyperfineA::hfsA(hw.get(), phi);
    if (rpa) {
      auto a = DiracOperator::HyperfineA::convertRMEtoA(phi, phi);
      Ap += a * rpap->dV(phi, phi);
      Ab += a * rpab->dV(phi, phi);
      Aw += a * rpaw->dV(phi, phi);
    }
    auto Fball = ((Ab / Ap) - 1.0) * 100.0;
    auto Fbw = ((Aw / Ap) - 1.0) * 100.0;
    bw.emplace_back(phi.n, phi.k, Fball, Fbw);
    printf("%7s| %12.5e %12.5e %12.5e | %9.5f  %9.5f \n", phi.symbol().c_str(),
           Ap, Ab, Aw, Fball, Fbw);
  }

  const auto do_screening = input.get("screening", false);
  if (do_screening) {
    // Create H-like orbitals: note: uses same grid
    std::vector<ScreenBW> H_bw;
    auto Hlike = Wavefunction(wf.rgrid->params(), wf.get_nuclearParameters(),
                              wf.alpha / PhysConst::alpha);
    std::cout << "\nH-like:\n" << Hlike.rgrid->gridParameters() << "\n";
    Hlike.localValence(DiracSpinor::state_config(wf.valence));
    Hlike.printValence();

    std::cout << "\nTabulate A (Mhz), and Bohr-Weisskopf effect eps(%): "
              << Hlike.atom() << " [H-like]"
              << "\n       |A:      Point         Ball           SP |e:    "
                 "Ball         SP\n";
    for (const auto &phi : Hlike.valence) {
      auto Ap = HyperfineA::hfsA(hp.get(), phi);
      auto Ab = HyperfineA::hfsA(hb.get(), phi);
      auto Aw = HyperfineA::hfsA(hw.get(), phi);
      auto Fball = ((Ab / Ap) - 1.0) * 100.0;
      auto Fbw = ((Aw / Ap) - 1.0) * 100.0;
      H_bw.emplace_back(phi.n, phi.k, Fball, Fbw);
      printf("%7s| %12.5e %12.5e %12.5e | %9.5f  %9.5f \n",
             phi.symbol().c_str(), Ap, Ab, Aw, Fball, Fbw);
    }

    std::cout << "\nScreening factors (ball, SP)\n";
    for (const auto &n : bw) {
      std::cout << n.n << " " << n.k << ":\n";
      for (const auto &m : H_bw) {
        if (m.k != n.k)
          continue;
        printf("    %2i %.5f %.5f\n", m.n, n.b / m.b, n.sp / m.sp);
      }
    }
  }
}

//******************************************************************************
void BW_eta_sp(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;

  input.checkBlock_old({});

  auto options = input.getBlock("options");
  auto sub_input = IO::InputBlock("hfs", {});
  if (options) {
    sub_input.add(options->options());
  }
  // sub_input.print();

  std::cout << "\nScreening: eta_sp = eps{ns}/eps{(n+1)p}\n"
            << "\nNo RPA (mainly for H-like!)\n";

  auto point_in = sub_input;
  auto ball_in = sub_input;
  auto BW_in = sub_input;

  point_in.add("F(r)=pointlike;");
  ball_in.add("F(r)=ball;");
  const auto doubly_odd = (wf.Anuc() % 2 == 0);
  if (doubly_odd)
    BW_in.add("F(r)=doublyOddBW;");
  else
    BW_in.add("F(r)=VolotkaBW;");

  auto hp = generateOperator(point_in, wf, false);
  auto hb = generateOperator(ball_in, wf, false);
  auto hw = generateOperator(BW_in, wf);

  std::cout << "\n      A0(MHz)         e(ball)   e(SP)     eta(b)   eta(sp)\n";
  for (const auto &Fs : wf.valence) {
    if (Fs.k != -1)
      continue;
    auto Asp = HyperfineA::hfsA(hp.get(), Fs);
    auto Asb = HyperfineA::hfsA(hb.get(), Fs);
    auto Asw = HyperfineA::hfsA(hw.get(), Fs);

    auto es_b = 100.0 * (Asb - Asp) / Asp;
    auto es_w = 100.0 * (Asw - Asp) / Asp;

    printf("%4s  %.7e  %7.5f  %7.5f\n", Fs.shortSymbol().c_str(), Asp, es_b,
           es_w);

    auto Fp = wf.getState(Fs.n + 1, +1);
    if (Fp) {

      auto App = HyperfineA::hfsA(hp.get(), *Fp);
      auto Apb = HyperfineA::hfsA(hb.get(), *Fp);
      auto Apw = HyperfineA::hfsA(hw.get(), *Fp);

      auto ep_b = 100.0 * (Apb - App) / App;
      auto ep_w = 100.0 * (Apw - App) / App;

      auto eta_sp_b = es_b / ep_b;
      auto eta_sp_w = es_w / ep_w;

      printf("%4s  %.7e  %7.5f  %7.5f  %7.4f  %7.4f\n",
             Fp->shortSymbol().c_str(), App, ep_b, ep_w, eta_sp_b, eta_sp_w);
    }
  }
}

} // namespace Module
