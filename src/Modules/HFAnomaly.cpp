#include "Modules/HFAnomaly.hpp"
#include "DiracOperator/Operators.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "MBPT/DiagramRPA.hpp"
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
void HFAnomaly(const IO::UserInputBlock &input, const Wavefunction &wf) {

  input.checkBlock({"rpa", "mu", "I", "rrms", "parity", "l", "gl", "mu1", "gl1",
                    "l1", "l2", "I1", "I2", "A"});

  const auto rpa = input.get("rpa", false);
  const auto Alist = input.get_list("A", std::vector<int>{});

  const auto sub_input = input.subBlock("MatrixElements::hfs",
                                        {"mu", "I", "rrms", "parity", "l", "gl",
                                         "mu1", "gl1", "l1", "l2", "I1", "I2"});
  const auto point_in = sub_input.copy_with("F(r)=pointlike");

  const auto ball_in = sub_input.copy_with("F(r)=ball");
  const auto doubly_odd = (wf.Anuc() % 2 == 0);
  const auto sp_in = doubly_odd ? sub_input.copy_with("F(r)=doublyOddBW")
                                : sub_input.copy_with("F(r)=VolotkaBW");

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
  const auto mu0 = input.get("mu", isotope0.mu);
  const auto I0 = input.get("I", isotope0.I_N);
  const auto r_rms0 = input.get("rrms", isotope0.r_rms);
  const auto gI0 = mu0 / I0;

  // RPA:
  std::unique_ptr<MBPT::DiagramRPA> rpap{nullptr}, rpab{nullptr}, rpas{nullptr};
  if (rpa) {
    std::cout << "\nIncluding RPA (diagram method) - must have basis\n";
    rpap = std::make_unique<MBPT::DiagramRPA>(hpt.get(), wf.basis, wf.core,
                                              wf.identity());
    rpab = std::make_unique<MBPT::DiagramRPA>(hbl.get(), rpap.get());
    rpas = std::make_unique<MBPT::DiagramRPA>(hsp.get(), rpap.get());
    std::cout << "Solving RPA core for point, ball, SP:\n";
    rpap->rpa_core(0.0, true);
    rpab->rpa_core(0.0, true);
    rpas->rpa_core(0.0, true);
  }

  // Calculate HFS for each state of reference isotope: including BW effect
  std::vector<ThreeAs> As;
  std::cout << "\nReference: A = " << wf.Anuc() << ", gI = " << gI0
            << ", r_rms = " << r_rms0 << "\n\n";
  std::cout
      << "      A(point)     A(ball)      A(SP)       | e(ball)  e(SP) [%]\n";
  for (const auto &Fv : wf.valence) {
    auto point = DiracOperator::Hyperfine::hfsA(hpt.get(), Fv);
    auto ball = DiracOperator::Hyperfine::hfsA(hbl.get(), Fv);
    auto sp = DiracOperator::Hyperfine::hfsA(hsp.get(), Fv);
    if (rpa) {
      const auto a = DiracOperator::Hyperfine::convertRMEtoA(Fv, Fv);
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
    wfA.hartreeFockCore("HartreeFock", 0.0, wf.coreConfiguration_nice());
    wfA.hartreeFockValence(DiracSpinor::state_config(wf.valence));
    wfA.basis = wf.basis; // OK??

    const auto hpt2 = generateOperator(
        {"MatrixElements::hfs", "F(r)=pointlike;"}, wfA, false);
    const auto hbl2 =
        generateOperator({"MatrixElements::hfs", "F(r)=ball;"}, wfA, false);
    const auto hsp2 = generateOperator(
        {"MatrixElements::hfs", "F(r)=VolotkaBW;"}, wfA, false);

    std::unique_ptr<MBPT::DiagramRPA> rpap2{nullptr}, rpab2{nullptr},
        rpas2{nullptr};
    if (rpa) {
      std::cout << "Including RPA (diagram method) - must have basis\n";
      // OK ? or need run with new basis!?
      rpap2 = std::make_unique<MBPT::DiagramRPA>(hpt2.get(), wfA.basis,
                                                 wfA.core, wfA.identity());
      rpab2 = std::make_unique<MBPT::DiagramRPA>(hbl2.get(), rpap.get());
      rpas2 = std::make_unique<MBPT::DiagramRPA>(hsp2.get(), rpap.get());
      std::cout << "Solving RPA core for point, ball, SP:\n";
      rpap2->rpa_core(0.0, true);
      rpab2->rpa_core(0.0, true);
      rpas2->rpa_core(0.0, true);
    }

    // for (const auto &Fv : wf.valence)
    std::cout << "\nA = " << nuc.A << ", gI = " << gI
              << ", r_rms = " << nuc.r_rms << "\n";
    std::cout << "      A(point)      e(ball)   e(SP)   | 1D2(ball) "
                 "1D2(SP) [%]\n";
    for (std::size_t i = 0; i < wfA.valence.size(); ++i) {
      const auto &Fv = wfA.valence[i];
      const auto [pt0, bl0, sp0] = As[i];
      auto point = DiracOperator::Hyperfine::hfsA(hpt2.get(), Fv);
      auto ball = DiracOperator::Hyperfine::hfsA(hbl2.get(), Fv);
      auto sp = DiracOperator::Hyperfine::hfsA(hsp2.get(), Fv);
      if (rpa) {
        const auto a = DiracOperator::Hyperfine::convertRMEtoA(Fv, Fv);
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
void HF_rmag(const IO::UserInputBlock &input, const Wavefunction &wf) {
  // For isotope 1 and 2
  // Loops over many values for magnetic radius of isotope 1, Rmag(1).
  // For each, finds Rmag(2) that reproduces a given hyperfine anomaly.
  // Potential to use two states to find rmag. If not, still interesting?

  std::cout << "\nTuning Rmag to fit hyperfine anomaly\n";

  input.checkBlock({"n", "kappa", "A2", "1D2", "rpa", "num_steps", "mu1", "mu2",
                    "I1", "I2", "eps_targ", "e1", "e2"});

  // A(1) is wf
  // A(2) is wf2
  const auto A2 = input.get<int>("A2");
  const auto eps_t = input.get("eps_targ", 1.0e-4);

  const auto n = input.get<int>("n");
  const auto kappa = input.get("kappa", -1);

  auto wf2 = Wavefunction(wf.rgrid->params(), {wf.Znuc(), A2},
                          wf.alpha / PhysConst::alpha);
  wf2.hartreeFockCore("HartreeFock", 0.0, wf.coreConfiguration_nice());
  wf2.hartreeFockValence(DiracSpinor::state_config(wf.valence));
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

  std::unique_ptr<MBPT::DiagramRPA> rpa01{nullptr}, rpa02{nullptr},
      rpa1{nullptr}, rpa2{nullptr}, rpa2a{nullptr}, rpa2b{nullptr};

  const double x = 0.1;
  const auto num_steps = input.get("num_steps", 40);
  const auto dr = 3.0 * x * rN0_au / (num_steps + 1);

  const auto h01 = DiracOperator::Hyperfine(
      mu1, I1, 0.0, *wf.rgrid, DiracOperator::Hyperfine::pointlike_F());
  const auto h02 = DiracOperator::Hyperfine(
      mu2, I2, 0.0, *wf.rgrid, DiracOperator::Hyperfine::pointlike_F());
  double dv01 = 0.0;
  double dv02 = 0.0;
  if (rpa) {
    rpa01 = std::make_unique<MBPT::DiagramRPA>(&h01, wf.basis, wf.core,
                                               wf.identity());
    rpa02 = std::make_unique<MBPT::DiagramRPA>(&h02, wf.basis, wf.core,
                                               wf.identity());
    rpa01->rpa_core(0.0);
    rpa02->rpa_core(0.0);
    auto a = DiracOperator::Hyperfine::convertRMEtoA(F1v, F1v);
    dv01 = a * rpa01->dV(F1v, F1v);
    dv02 = a * rpa02->dV(F1v, F1v);
  }
  double a0_1 = h01.hfsA(F1v) + dv01;
  double a0_2 = h02.hfsA(F1v) + dv02;

  const auto d12_targ = input.get<double>("1D2", 0.0);
  if (d12_targ != 0.0) {
    std::cout << "\nA0(1) = " << a0_1 << "\n";
    std::cout << "A0(2) = " << a0_2 << "\n";
    std::cout << "1D2 target: " << d12_targ << "\n";
    std::cout << "R0(1) = " << rN0_au * PhysConst::aB_fm << "\n";
    std::cout << "R0(2) = " << rN2 * PhysConst::aB_fm << "\n";
    std::cout << "Rmag(1)  Rmag(2)  e1     e2     1D2    eps(D)   del(R)\n";
    for (double rN = (1.0 - x) * rN0_au; rN < (1.0 + 2 * x) * rN0_au;
         rN += dr) {

      const auto h1 = DiracOperator::Hyperfine(mu1, I1, rN, *wf.rgrid, Fr1);
      double dv1 = 0.0;
      if (rpa) {
        rpa1 = std::make_unique<MBPT::DiagramRPA>(&h1, rpa01.get());
        rpa1->grab_tam(rpa01.get()); // don't start from scratch
        rpa1->rpa_core(0.0, false);
        auto a = DiracOperator::Hyperfine::convertRMEtoA(F1v, F1v);
        dv1 = a * rpa1->dV(F1v, F1v);
      }

      const auto a1 = h1.hfsA(F1v) + dv1;

      // Use halving-interval method to find rmag(2) that reproduces 1D2, given
      // rmag(1)
      auto r2a = (1.0 - x) * rN2;
      auto r2 = rN2;
      auto r2b = (1.0 + x) * rN2;
      int tries = 0;
      double eps;
      double del_r;
      double d12;
      double bw1, bw2;
      while (tries++ < 40) {

        // Note: For RPA, could make this 3x faster!
        // Only need to change the 'outer' guess each time, do not need to
        // re-evaluate RPA for the other two guesses! Complicated though??
        const auto h2a = DiracOperator::Hyperfine(mu2, I2, r2a, *wf.rgrid, Fr2);
        const auto h2 = DiracOperator::Hyperfine(mu2, I2, r2, *wf.rgrid, Fr2);
        const auto h2b = DiracOperator::Hyperfine(mu2, I2, r2b, *wf.rgrid, Fr2);

        double dv2 = 0.0;
        double dv2a = 0.0;
        double dv2b = 0.0;
        if (rpa) {
          rpa2a = std::make_unique<MBPT::DiagramRPA>(&h2a, rpa1.get());
          rpa2 = std::make_unique<MBPT::DiagramRPA>(&h2, rpa1.get());
          rpa2b = std::make_unique<MBPT::DiagramRPA>(&h2b, rpa1.get());
          rpa2->grab_tam(rpa02.get());  // don't start from scratch
          rpa2a->grab_tam(rpa02.get()); // don't start from scratch
          rpa2b->grab_tam(rpa02.get()); // don't start from scratch
          rpa2a->rpa_core(0.0, false);
          rpa2->rpa_core(0.0, false);
          rpa2b->rpa_core(0.0, false);
          const auto a = DiracOperator::Hyperfine::convertRMEtoA(F2v, F2v);
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

        // std::cout << tries << " " << r2a * PhysConst::aB_fm << " "
        //           << r2b * PhysConst::aB_fm << " " << d12 << "\n";

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
      printf("%.5f  %.5f %6.3f %6.3f %7.4f %.1e  %.1e  [%i]\n",
             rN * PhysConst::aB_fm, r2 * PhysConst::aB_fm, bw1, bw2, d12, eps,
             del_r * PhysConst::aB_fm, tries);
      // std::cout << bw1 << " " << bw2 << "\n";
    } // rMag
  }

  std::cout << "\n";
  auto etarg1 = input.get("e1", 0.0);
  if (etarg1 != 0.0) {
    std::cout << "Fitting rmag to find epsilon for A1: e_targ=" << etarg1
              << "\n";
    calc_thing(F1v, etarg1, rN0_au, mu1, I1, l1, gl);
  }
  std::cout << "\n";
  auto etarg2 = input.get("e2", 0.0);
  if (etarg2 != 0.0) {
    std::cout << "Fitting rmag to find epsilon for A2: e_targ=" << etarg2
              << "\n";
    calc_thing(F2v, etarg2, rN2, mu2, I2, l2, gl);
  }
}

//******************************************************************************
static void calc_thing(const DiracSpinor &Fv, double e_targ, double r0,
                       double mu, double I, int l, int gl) {
  const auto Fr = DiracOperator::Hyperfine::volotkaBW_F(mu, I, l, gl);

  const auto h0 = DiracOperator::Hyperfine(
      mu, I, 0.0, *Fv.rgrid, DiracOperator::Hyperfine::pointlike_F());
  double dv0 = 0.0; // XXX RPA!
  const auto A0 = h0.hfsA(Fv) + dv0;

  std::cout << "nb: no RPA\n";

  int tries = 0.0;
  double r = r0;
  double ra = 0.8 * r0, rb = 1.2 * r0;
  while (tries++ < 100) {
    const auto ha = DiracOperator::Hyperfine(mu, I, ra, *Fv.rgrid, Fr);
    const auto h = DiracOperator::Hyperfine(mu, I, r, *Fv.rgrid, Fr);
    const auto hb = DiracOperator::Hyperfine(mu, I, rb, *Fv.rgrid, Fr);

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
      auto c = PhysConst::aB_fm;
      // std::cout << r * c << " | " << bw << " eps_r=" << (rb - ra) / r
      //           << " eps_e=" << (bwa - bwb) / bw << "\n";
      printf("rN=%8.6f, e=%6.3f,   eps_r=%.0e  eps_e=%.0e\n", r * c, bw,
             (rb - ra) / r, (bwa - bwb) / bw);
      break;
    }
  }
}

//******************************************************************************
void calculateBohrWeisskopf(const IO::UserInputBlock &input,
                            const Wavefunction &wf) {
  using namespace DiracOperator;

  input.checkBlock({"rpa", "rpa_diagram", "screening", "mu", "I", "rrms",
                    "F(r)", "parity", "l", "gl", "mu1", "gl1", "l1", "l2", "I1",
                    "I2", "printF"});

  IO::UserInputBlock point_in("MatrixElements::hfs", input);
  IO::UserInputBlock ball_in("MatrixElements::hfs", input);
  IO::UserInputBlock BW_in("MatrixElements::hfs", input);
  point_in.add("F(r)=pointlike");
  ball_in.add("F(r)=ball");
  if (wf.Anuc() % 2 == 0)
    BW_in.add("F(r)=doublyOddBW");
  else
    BW_in.add("F(r)=VolotkaBW");

  auto hp = generateOperator(point_in, wf, false);
  auto hb = generateOperator(ball_in, wf, false);
  auto hw = generateOperator(BW_in, wf);

  // nb: can only do diagram RPA for hfs
  const auto rpa = input.get("rpa", false) || input.get("rpa_diagram", false);

  std::unique_ptr<MBPT::DiagramRPA> rpap{nullptr}, rpab{nullptr}, rpaw{nullptr};
  if (rpa) {
    std::cout << "\nIncluding RPA (diagram method) - must have basis\n";
    rpap = std::make_unique<MBPT::DiagramRPA>(hp.get(), wf.basis, wf.core,
                                              wf.identity());
    rpab = std::make_unique<MBPT::DiagramRPA>(hb.get(), rpap.get());
    rpaw = std::make_unique<MBPT::DiagramRPA>(hw.get(), rpap.get());
    std::cout << "Solving RPA core for point, ball, SP:\n";
    rpap->rpa_core(0.0);
    rpab->rpa_core(0.0);
    rpaw->rpa_core(0.0);
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
    auto Ap = Hyperfine::hfsA(hp.get(), phi);
    auto Ab = Hyperfine::hfsA(hb.get(), phi);
    auto Aw = Hyperfine::hfsA(hw.get(), phi);
    if (rpa) {
      auto a = DiracOperator::Hyperfine::convertRMEtoA(phi, phi);
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
      auto Ap = Hyperfine::hfsA(hp.get(), phi);
      auto Ab = Hyperfine::hfsA(hb.get(), phi);
      auto Aw = Hyperfine::hfsA(hw.get(), phi);
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

} // namespace Module
